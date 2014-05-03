from subprocess import PIPE, Popen, CalledProcessError
import re
import json
import time
import datetime

SAMPLE_PROCESS_PROGRAMS = ["pairs", "seqqs", "scythe", "seqtk", "bwa", "samtools"]

SAMPLE_PROCESS_CMD = """\
set -o pipefail; \
{pairs} join {read1} {read2} | \
{seqqs} -i -p {stats_dir}/raw_{sample_id} -e - | \
{scythe} -a {adapters_file} -p {prior} - | \
{seqtk} trimfq -q {error} - | \
{seqqs} -i -p {stats_dir}/processed_{sample_id} -e - | \
{bwa} mem -M -t {nthreads} -R '@RG\tID:{read_id}\tPL:illumina\tSM:{sample_id}' -v 1 -p {reference} - | \
{samtools} view -b -S -u - | samtools sort -@ {nthreads} -m {mem} - {bam_dir}/{prefex}.sorted
"""

SLURM_BATCH = """\
#SBATCH -o {log_dir}/{jobname}-%%j-%%a.stdout
#SBATCH -e {log_dir}/{jobname}-%%j-%%a.stderr
#SBATCH -J {jobname}
#SBATCH --cpus-per-task={nthreads}
#SBATCH --mem-per-cpu={mem}
#SBATCH --array=1-{nsamples}

module load bwa samtools seqqs scythe

sed -n "$SLURM_ARRAY_TASK_ID"p {sample_config} | python {sample_dispatch_py} runner
"""

def get_template_keys(cmd):
    return re.findall(r'{(\w+)}', cmd)

def safe_templater(cmd, mapping):
    """
    A safer version of str.format, but checks that all keys are being used.
    """
    template_keys = get_template_keys(cmd)
    key_symdiff = set(template_keys).symmetric_difference(set(mapping))
    if length(key_symdiff):
        raise ValueError("command template's keys are not identical to keys in mapping: " + ', '.join(key_symdiff))
    return cmd.format(**mapping)

def validate_setupfile(setup_params):
    """
    Validate setup by checking that programs can be found.
    """
    for command, path in setup_params.items():
        try:
            check_call("command -v %s" % path)
        except CalledProcessError:
            raise CalledProcessError("command '%s' with path '%s' not found." % (cmd, path))

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

def validate_directory(direc):
    if not os.path.isdir(direc):
        raise IOError("directory '%s' does not exist." % direct)

def dispatch(arg):
    """
    Create a sample JSON file for a run, and launch it using Popen.
    """
    # validate that programs in setup file exist
    setup_params = json.load(args.setup)
    validate_setupfile(setup_params)

    # validate reference
    validate_reference(args.ref)

    # validate directories
    validate_directory(args.log)
    validate_directory(args.stats)

    # create batch script
    sbatch_params = {"log_dir":args.log, "jobname":args.job, "nthreads":nthreads,
                    "mem":args.mem, "nsamples":nsamples, "sample_dspatch_py":__file__,
                    "sample_config":args.run}
    batch_script = safe_templater(SLURM_BATCH, sbatch_params)
    batch_file = "%s_batch.sh" % args.run
    with open(batch_file, 'w') as f:
        f.write(batch_script)

    # create sample config JSON file, starting off with global config passed through args
    global_sample_config = dict(ref=arg.ref, adapter=args.adapter, prior=str(args.prior), error=args.trim_error,
                         log_dir=args.log, stats_dir=args.stats, reference=args.ref, nthreads=args.threads,
                         mem=args.mem).items()
    sample_config = dict()
    config_file = open("%s_samples.txt", 'w')
    with open(args.samples) as f:
        for line in f:
            sample_id, read_id, reads1, reads2 = line.strip().split("\t")
            this_sample_config = dict(sample_id=sample_id, read_id=read_id, reads1=reads1, reads2=reads2)
            combined_sample_config = dict(global_sample_config.items() + this_sample_config.items())
            sample_config[sample_id] = combined_sample_config
            config_file.write(json.dumps(combined_sample_config))

    if not args.dry_run:
        # now, start the batch script
        p = Popen(["sbatch", batch_file])
        p.wait() # this should be fast.
        if p.returncode != 0:
            sys.exit("step dispatch exited abnormally with return code %d" % (p.returncode))

def runner(args):
    """
    Run a sample through an NGS pipeline (as a command) using Popen.
    TOOD: logging, timing.
    """
    sample_params = json.loads(sys.stdin.readline().rstrip())
    setup_params = json.load(args.setup)
    cmd = sample_templater(SAMPLE_PROCESS_CMD, sample_params)
    sys.stderr.write("[%s] starting preprocessing and alignment of sample '%s'" %
                      (datetime.now().isoformat(' '), sample_params['sample_id']))
    if args.dry_run:
        sys.stdout.write("[debug] command for '%s': %s" % (sample_params['sample_id'], cmd))
        return
    tstart = time.time()
    p = Popen(cmd)
    p.wait()
    if p.returncode != 0:
        # make this as loud as possible so Slurm can handle it
        sys.exit("sample '%s' exited abnormally with return code %d" %
                (setup_params['sample_id'], p.returncode))
    tend = time.time()
    elapsed = tend - tstart
    sys.stderr.write("[%s] completed preprocessing and alignment of sample '%s' in %s seconds" %
                      (datetime.now().isoformat(' '), sample_params['sample_id'], str(round(elapsed, 5))))

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Dispatch a set of samples to be preprocessed and aligned.')
    parser.add_argument('--dry-run', '-d', help="don't run batch file", action="store_true")
    subparsers = parser.add_subparsers(help='sub-command help')
    parser_dispatch = subparsers.add_parser('dispatch', help='create a set of samples in JSON format to dispatch')
    parser_dispatch.add_argument('--ref', '-r', help="reference FASTA file, indexed by BWA")
    parser_dispatch.add_argument('--samples', '-s', help="tab-delimited sample sheet")
    parser_dispatch.add_argument('--adapter', '-a', help="adapter file (for scythe")
    parser_dispatch.add_argument('--prior', '-p', help="prior adapter contamination rate (for scythe)",
                                type=float, default=0.3)
    parser_dispatch.add_argument('--trim-error', '-e', help="trimmer error rate threshold (for seqtk trimfq)",
                                type=float, default=0.05)
    parser_dispatch.add_argument('--log', '-l', help="directory for logging", default="log/")
    parser_dispatch.add_argument('--stats', '-S', help="directory for diagnostic statistics files", default="stats/")
    parser_dispatch.add_argument('--mem', '-m', help="memory (in MB) to use *per* thread in sorting", type=int)
    parser_dispatch.add_argument('--threads', '-t', help="threads to use in alignment and sorting", default=1, type=int)
    parser_dispatch.add_argument('--job', '-j', help="job name")
    parser_dispatch.add_argument('run', help="name of a run")
    parser_dispatch.set_defaults(func=dispatch)

    parser_runner = subparsers.add_parser('runner', help="take a single sample JSON entry from standard in and process")
    parser_runner.add_argument('--setup', '-S', help="setup JSON file (with file paths to programs)")
    parser_runner.set_defaults(func=runner)

    args = parser_dispatch.parse_args()
    args.func(args)


