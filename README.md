# SNP Pipeline

This set of scripts implements a simple read preprocessing and alignment
pipeline using Slurm array jobs. Many steps are configurable, but the general
pipeline is:

1. Join reads into interleaved pairs (to simplify processing).

2. Run reads through [seqqs](https://github.com/vsbuffalo/seqqs), which records metrics on read quality, length, base composition.

3. Trim adapter sequences off of reads using [scythe](https://github.com/vsbuffalo/scythe).

4. Quality-based trimming with [seqtk's](https://github.com/lh3/seqtk) trimfq.

5. Another around of seqqs, which records post pre-processing read quality metrics.

6. Align with [BWA-MEM](https://github.com/lh3/bwa). The exact command is:

        mem -M -t <nthreads> -R <readgroup> -v 1 -p <reference> <in.fq>

7. Convert reads to BAM and sort reads by position with [samtools](https://github.com/samtools/samtools).

This pipeline is just a template; it can be adjusted easily in code. In the
future, this template may be factored out of code.

The pipeline in two steps (which are hidden from the user):

First, the *dispatch* step creates a text file of all the sample's run
information (called the `sample_config.json` file). Each of these is a JSON
entry in a text file, which is passed to a runner. Each of these JSON entries
contains all data needed by the processing command: reference, parameters,
sample information, etc.

Then, the pipeline writes a Slurm batch script in your directory for this job,
and runs it. This batch script then calls the *runners* (using Slurm's array
job infrastructure) which is a subcommand of this pipeline that takes a single
line from the `sample_config.json` file and run it.

## Running the SNP Pipeline

All processing is done on a single pair of reads. In order to run, two
configuration files are needed (note you can name these whatever you want):

1. A `setup.json` JSON configuration file, which includes paths to all programs.

2. A `samples.txt` file which provides a mapping between both read pairs,
	 sample IDs (`@RG:SM` in the SAM spec), and read group IDs (`@RG:ID` in the
   SAM spec). The header format is:

        sample_id, read_id, reads1, reads2





