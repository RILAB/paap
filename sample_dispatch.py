import os
import argparse
import sys
from subprocess import PIPE, Popen, CalledProcessError, check_call, check_output
import re
import json
import datetime
from logbook import Logger
from collections import OrderedDict

DEVNULL = open(os.devnull, 'wb')

### Pipeline Components ###

## Pipeline Initilization
PIPE_INITIALIZATION_STEPS = OrderedDict([
        ("pipefail", "set -o pipefail; "),
        ("interleave", "{pairs} join {reads1} {reads2} | ")])

## Pre-processing Steps
# these are optional; keys correspond to CL args
PREPROCESSING_STEPS = OrderedDict([
        ("pre_seqqs", "{seqqs} -i -p {stats_dir}/raw_{sample_id} -e - | "),
        ("scythe", "{scythe} -a {adapters_file} -p {prior} - 2> {stats_dir}/scythe_{sample_id}.txt | "),
        ("trimfq", "{seqtk} trimfq -q {error} - | "),
        ("post_seqqs", "{seqqs} -i -p {stats_dir}/processed_{sample_id} -e - | ")])

## Alignment Steps and Post-Alignment Steps
# alignment and post_alignment steps required
ALN_STEPS = OrderedDict([
        ("bwa", "{bwa} mem -M -t {nthreads} -R '@RG\tID:{read_id}\tPL:illumina\tSM:{sample_id}' -v 1 -p {reference} - | ")])
POST_ALN_STEPS = OrderedDict([
        ("samtools-to-bam", "{samtools} view -b -S -u - > {bam_dir}/{sample_id}.bam")])

## Sort step (run not in pipe, as this causes memory issues)
SORT_STEPS = OrderedDict([
        ("samtools-sort", "samtools sort -@ {nthreads} -m {mem}M {bam_dir}/{sample_id}.bam {bam_dir}/{sample_id}_sorted")])

SLURM_BATCH = """\
#!/bin/bash
#SBATCH -o {log_dir}/{jobname}-%j-%a.stdout
#SBATCH -e {log_dir}/{jobname}-%j-%a.stderr
#SBATCH -J {jobname}
#SBATCH --cpus-per-task={nthreads}
#SBATCH --mem-per-cpu={mem}
#SBATCH --array=1-{nsamples}
#SBATCH --partition={partition}

module load bwa samtools seqqs scythe

sed -n "$SLURM_ARRAY_TASK_ID"p {sample_config} | python {sample_dispatch_py} runner
"""

def merge_steps(collection, steps):
    """
    Given a collection of processing steps, merge those specified
    by steps, creating a command. No checking; it's responsibility
    of collection to ensure parts are interoperable.
    """
    steps = set(steps)
    parts = [step.strip() for key, step in collection.items() if key in steps]
    return " ".join(parts)

def build_sample_aln_command(sample_params):
    """
    Construct a read pair preprocessing and alignment command from
    templates and ordered steps.
    """
    steps = list()
    steps.append(merge_steps(PIPE_INITIALIZATION_STEPS,
                            ("pipefail", "interleave")))

    # preprocessing includes many optional components
    steps.append(merge_steps(PREPROCESSING_STEPS, sample_params["preprocess-steps"]))

    # these components are (so far) not optional
    steps.append(ALN_STEPS["bwa"])
    steps.append(POST_ALN_STEPS["samtools-to-bam"])
    cmd_template = " ".join(steps)
    return safe_templater(cmd_template, sample_params)

def validate_program_exists(command):
    p = check_call("command -v %s" % command, shell=True, stdout=DEVNULL, stderr=DEVNULL)
    if p != 0:
        raise ValueError("program '%s' not found." % command)

def find_bash():
    try:
        bash_path = check_output(["which", "bash"])
    except CalledProcessError:
        raise ValueError("cannot find bash")
    validate_program_exists(bash_path)
    return bash_path.rstrip()

def get_template_keys(cmd):
    return re.findall(r'{(\w+)}', cmd)

def safe_templater(cmd, mapping):
    """
    A safer version of str.format, but checks that no template keys are unfilled.
    """
    template_keys = get_template_keys(cmd)
    key_diff = set(template_keys) - set(mapping)
    if len(key_diff):
        raise ValueError("command template's keys contain keys in mapping: " + ', '.join(key_diff))
    return cmd.format(**mapping)

def validate_setupfile(setup_params):
    """
    Validate setup by checking that programs can be found.
    """
    for command, path in setup_params.items():
        try:
            validate_program_exists(command)
        except:
            raise ValueError("command '%s' with path '%s' not found." % (command, path))

def validate_reference(reference):
    """
    Ensure reference exists and has the necessary BWA index files by comparing
    expected suffixes.
    """
    bwa_index_suffixes = ".amb .ann .bwt .pac .sa".split()
    expected_index_files = [reference] + ["%s%s" % (reference, suf) for suf in bwa_index_suffixes]
    for index_file in expected_index_files:
        if not os.path.isfile(index_file):
            raise IOError("expected file '%s' required for BWA alignment does not exist." % index_file)

def validate_directory(direc, logger):
    if not os.path.isdir(direc):
        logger.info("directory '%s' does not exist, creating it." % direc)
        os.makedirs(direc)

def create_sample_config(sample_file, json_out_file, global_params):
    """
    Parse the tab-delimited sample file, saving each run configuration as a
    JSON entry in `json_out_file`. Return a dictionary of all sample
    configurations, built off `global_params`.
    """
    sample_config = dict()
    config_file = open(json_out_file, 'w')
    for line in sample_file:
        sample_id, read_id, reads1, reads2 = line.strip().split("\t")
        reads_exist = dict([(read, os.path.isfile(read)) for read in (reads1, reads2)])
        if not all(reads_exist.values()):
            msg = ', '.join(["'%s'" % read for read, found in reads_exist.items() if not found])
            logger.critical("reads file(s) %s not found." % msg)
            sys.exit(1)
        this_sample_config = dict(sample_id=sample_id, read_id=read_id, reads1=reads1, reads2=reads2)
        combined_sample_config = dict(global_params.items() + this_sample_config.items())
        sample_config[sample_id] = combined_sample_config
        config_file.write(json.dumps(combined_sample_config) + "\n")
    return sample_config

def dispatch(args):
    """
    Create a sample JSON file for a run, and launch it using Popen.
    """
    dispatch_log = Logger('dispatch')
    # validate that programs in setup file exist
    setup_params = json.load(args.setup)
    validate_setupfile(setup_params)

    # validate reference
    validate_reference(args.ref)

    # validate directories
    validate_directory(args.log, dispatch_log)
    validate_directory(args.stats, dispatch_log)
    validate_directory(args.bam_dir, dispatch_log)

    if (args.scythe and args.adapter is None) or (args.adapter is not None and not os.path.isfile(args.adapter)):
        logger.critical("adapter file for Scythe no specified, or does not exist.")
        sys.exit(1)

    # create sample config JSON file, starting off with global config passed through args
    global_sample_config = dict(reference=args.ref, adapters_file=args.adapter,
                                prior=str(args.prior), error=args.trim_error,
                                stats_dir=args.stats, nthreads=args.threads,
                                mem=args.mem, bam_dir=args.bam_dir)

    # which preprocess steps to use
    global_sample_config["preprocess-steps"] = list()
    for step in PREPROCESSING_STEPS:
        if step in args and args.__getattribute__(step):
            global_sample_config["preprocess-steps"].append(step)

    global_params = dict(global_sample_config.items() + setup_params.items())
    sample_config = "%s_samples.txt" % args.job
    samples = create_sample_config(args.samples, sample_config, global_params)

    # create batch script
    sbatch_params = {"log_dir":args.log, "jobname":args.job, "nthreads":args.threads,
                    "mem":args.mem, "nsamples":len(samples), "sample_dispatch_py":__file__,
                    "sample_config":sample_config, "partition":args.partition}
    batch_script = safe_templater(SLURM_BATCH, sbatch_params)
    batch_file = "%s_batch.sh" % args.job
    with open(batch_file, 'w') as f:
        f.write(batch_script)

    if not args.dry_run:
        # now, start the batch script
        dispatch_log.info("submitting sbatch script '%s'." % batch_file)
        sbatch_cmd = ["sbatch"]

        if args.email is not None:
            sbatch_cmd.extend(["--mail-type", "ALL"])
            sbatch_cmd.extend(["--mail-user", args.email])
        sbatch_cmd.append(batch_file)
        retcode = check_call(sbatch_cmd)
        if retcode != 0:
            dispatch_log.critical("submitting batch script '%s' exited abnormally with return code %d." % (batch_file, retcode))
            sys.exit(retcode.returncode)
        dispatch_log.info("submitting sbatch script '%s' complete." % batch_file)

def run_command_on_sample(cmd, logger, sample, desc):
    logger.info("%s starting %s" % (sample, desc))
    logger.info("%s command: %s" % (sample, cmd))
    tstart = datetime.datetime.now()
    p = Popen(cmd, shell=True, executable=find_bash())
    p.wait()
    if p.returncode != 0:
        # make this as loud as possible so Slurm can handle it
        logger.critical("%s exited abnormally with return code %d." % (sample, p.returncode))
        sys.exit(p.returncode)
    tend = datetime.datetime.now()
    elapsed = tend - tstart
    logger.info("%s completed %s in: %s" % (sample, desc, str(elapsed)))

def runner(args):
    """
    Run a sample through an NGS pipeline (as a command) using Popen.
    TOOD: logging, timing.
    """
    if args.config is None:
        sample_params = json.loads(sys.stdin.readline().rstrip())
    else:
        # read the first line from the test/debug config file
        sample_params = json.loads(args.config.readline().rstrip())

    sample = sample_params["sample_id"]
    runner_log = Logger("%s logger" % sample)
    tstart = datetime.datetime.now() # total run time

    # preprocessing and alignment
    aln_cmd = build_sample_aln_command(sample_params)
    run_command_on_sample(aln_cmd, runner_log, sample, desc="preprocessing and alignment")
    if args.dry_run:
        return

    sort_cmd = safe_templater(SORT_STEPS['samtools-sort'], sample_params)
    run_command_on_sample(sort_cmd, runner_log, sample, desc="sorting BAM file")
    tend = datetime.datetime.now()
    elapsed = tend - tstart
    runner_log.info("%s all processing completed in: %s." % (sample, str(elapsed)))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Dispatch a set of samples to be preprocessed and aligned.')
    parser.add_argument('--dry-run', '-n', help="don't run", action="store_true", default=False)
    subparsers = parser.add_subparsers(help='sub-command help')
    parser_dispatch = subparsers.add_parser('dispatch', help='create a set of samples in JSON format to dispatch')
    parser_dispatch.add_argument('--ref', '-r', help="reference FASTA file, indexed by BWA", required=True)

    # optional preprocessing steps
    preprocess = parser_dispatch.add_argument_group("pre-process",
            "steps to run during pre-processing")
    preprocess.add_argument('--trimfq', '-T', help="quality trim sequences with seqtk's trimfq",
                           default=False, action="store_true")
    preprocess.add_argument('--scythe', '-Y', help="trim adapter sequences with Scythe",
                           default=False, action="store_true")
    preprocess.add_argument('--pre-seqqs', '-q', default=False, action="store_true",
                           help="run seqqs before pre-processing to record statistics about quality")
    preprocess.add_argument('--post-seqqs', '-Q', default=False, action="store_true",
                           help="run seqqs after pre-processing to record statistics about quality")
    preprocess.add_argument('--adapter', '-a', help="adapter file (for scythe", default=None)
    preprocess.add_argument('--prior', '-p', help="prior adapter contamination rate (for scythe)",
                           type=float, default=0.3)
    preprocess.add_argument('--trim-error', '-e', help="trimmer error rate threshold (for seqtk trimfq)")

    # general configurations
    general_options = parser_dispatch.add_argument_group("general",
            "general options in pre-processing and alignment")
    general_options.add_argument('--samples', '-s', help="tab-delimited sample sheet", required=True,
                                type=argparse.FileType('r'))
    general_options.add_argument('--log', '-l', help="directory for logging", default="log/")
    general_options.add_argument('--stats', '-d',
                                help="directory for diagnostic statistics files", default="stats/")
    general_options.add_argument('--bam-dir', '-b', help="directory for output BAM files", required=True)
    general_options.add_argument('--mem', '-m', help="memory (in MB) to use *per* thread in sorting",
                                type=int, default=768)
    general_options.add_argument('--threads', '-t', help="threads to use in alignment and sorting",
                                default=1, type=int)
    general_options.add_argument('--setup', '-S', help="setup JSON file (with file paths to programs)",
                                required=True,
                                type=argparse.FileType('r'))
    general_options.add_argument('--job', '-j', help="job name", required=True)
    general_options.add_argument('--partition', '-P', help="Slurm partition", required=True)
    general_options.add_argument('--email', '-E', help="your email address for event notification",
            default=None)

    parser_dispatch.set_defaults(func=dispatch)

    # Runner: for internal use
    parser_runner = subparsers.add_parser('runner',
        help="interal use; runs samples serialized as JSON through pipeline")
    parser_runner.add_argument('config',
        help="config file argument (for debugging, normally passed through standard in)",
        default=None, nargs='?', type=argparse.FileType('r'))
    parser_runner.set_defaults(func=runner)

    args = parser.parse_args()
    args.func(args)


