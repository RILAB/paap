import os

import sys
from subprocess import PIPE, Popen, CalledProcessError, check_call, check_output
import re
import json
import time
from datetime import datetime
from logbook import Logger

DEVNULL = open(os.devnull, 'wb')

SAMPLE_PROCESS_PROGRAMS = ["pairs", "seqqs", "scythe", "seqtk", "bwa", "samtools"]

SAMPLE_PROCESS_CMD = """\
set -o pipefail; \
{pairs} join {reads1} {reads2} | \
{seqqs} -i -p {stats_dir}/raw_{sample_id} -e - | \
{scythe} -a {adapters_file} -p {prior} - | \
{seqtk} trimfq -q {error} - | \
{seqqs} -i -p {stats_dir}/processed_{sample_id} -e - | \
{bwa} mem -M -t {nthreads} -R '@RG\tID:{read_id}\tPL:illumina\tSM:{sample_id}' -v 1 -p {reference} - | \
{samtools} view -b -S -u - | samtools sort -@ {nthreads} -m {mem}M - {bam_dir}/{sample_id}.sorted \
"""

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
    A safer version of str.format, but checks that all keys are being used.
    """
    template_keys = get_template_keys(cmd)
    key_symdiff = set(template_keys).symmetric_difference(set(mapping))
    if len(key_symdiff):
        raise ValueError("command template's keys are not identical to keys in mapping: " + ', '.join(key_symdiff))
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

    # create sample config JSON file, starting off with global config passed through args
    global_sample_config = dict(reference=args.ref, adapters_file=args.adapter, prior=str(args.prior), error=args.trim_error,
                                stats_dir=args.stats, nthreads=args.threads,
                                mem=args.mem, bam_dir=args.bam_dir)
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
        retcode = check_call(["sbatch", batch_file])
        if retcode != 0:
            dispatch_log.critical("submitting batch script '%s' exited abnormally with return code %d." % (batch_file, retcode))
            sys.exit(retcode.returncode)
        dispatch_log.critical("submitting sbatch script '%s' complete." % batch_file)

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

    cmd = safe_templater(SAMPLE_PROCESS_CMD, sample_params)
    runner_log.info("%s starting preprocessing and alignment of sample." % sample)
    if args.dry_run:
        runner_log.debug("%s command: %s" % (sample, cmd))
        return
    tstart = time.time()
    p = Popen(cmd, shell=True, executable=find_bash())
    p.wait()
    if p.returncode != 0:
        # make this as loud as possible so Slurm can handle it
        runner_log.critical("%s exited abnormally with return code %d." % (sample, p.returncode))
        sys.exit(p.returncode)
    tend = time.time()
    elapsed = tend - tstart
    runner_log.info("%s completed preprocessing and alignment in %s seconds." % (sample, str(round(elapsed, 5))))

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Dispatch a set of samples to be preprocessed and aligned.')
    parser.add_argument('--dry-run', '-n', help="don't run", action="store_true", default=False)
    subparsers = parser.add_subparsers(help='sub-command help')
    parser_dispatch = subparsers.add_parser('dispatch', help='create a set of samples in JSON format to dispatch')
    parser_dispatch.add_argument('--ref', '-r', help="reference FASTA file, indexed by BWA", required=True)
    parser_dispatch.add_argument('--samples', '-s', help="tab-delimited sample sheet", required=True,
                                type=argparse.FileType('r'))
    parser_dispatch.add_argument('--adapter', '-a', help="adapter file (for scythe", required=True)
    parser_dispatch.add_argument('--prior', '-p', help="prior adapter contamination rate (for scythe)",
                                type=float, default=0.3)
    parser_dispatch.add_argument('--trim-error', '-e', help="trimmer error rate threshold (for seqtk trimfq)",
                                type=float, default=0.05)
    parser_dispatch.add_argument('--log', '-l', help="directory for logging", default="log/")
    parser_dispatch.add_argument('--stats', '-d', help="directory for diagnostic statistics files", default="stats/")
    parser_dispatch.add_argument('--bam-dir', '-b', help="directory for output BAM files", required=True)
    parser_dispatch.add_argument('--mem', '-m', help="memory (in MB) to use *per* thread in sorting", type=int, default=768)
    parser_dispatch.add_argument('--threads', '-t', help="threads to use in alignment and sorting", default=1, type=int)
    parser_dispatch.add_argument('--setup', '-S', help="setup JSON file (with file paths to programs)", required=True,
                                type=argparse.FileType('r'))
    parser_dispatch.add_argument('--job', '-j', help="job name", required=True)
    parser_dispatch.add_argument('--partition', '-P', help="Slurm partition", required=True)
    parser_dispatch.set_defaults(func=dispatch)

    parser_runner = subparsers.add_parser('runner', help="take a single sample JSON entry from standard in and process")
    parser_runner.add_argument('config', help="config file argument (for debugging, normally passed through standard in)",
                               default=None, nargs='?', type=argparse.FileType('r'))
    parser_runner.set_defaults(func=runner)

    args = parser.parse_args()
    args.func(args)


