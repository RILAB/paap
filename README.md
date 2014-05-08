# Preprocessing and Alignment Pipeline (PAAP)

This scripts implements a simple read preprocessing, alignment, and SNP calling
pipeline using Slurm array jobs. It's designed to be just a single Python
script, `paap.py` that you copy into project directory (and place under version
control). This way, there's never any question what steps you ran on your data
and what version of the script was used.


`paap.py` is (sadly yet another custom) Slurm NGS task processing -- use the
much better [bcbio-nextgen](https://github.com/chapmanb/bcbio-nextgen) for more
serious processing of model-organism data.

Many steps are configurable, but the general pipeline is:

**Pipeline initialization** 

1. Join reads into interleaved pairs (to simplify processing).

**Pre-processing**

2. Run reads through [seqqs](https://github.com/vsbuffalo/seqqs), which records metrics on read quality, length, base composition.

3. Trim adapter sequences off of reads using [scythe](https://github.com/vsbuffalo/scythe).

4. Quality-based trimming with [seqtk's](https://github.com/lh3/seqtk) trimfq.

5. Another around of seqqs, which records post pre-processing read quality metrics.

**Alignment**

6. Align with [BWA-MEM](https://github.com/lh3/bwa). The exact command is:

        bwa mem -M -t <nthreads> -R <readgroup> -v 1 -p <reference> <in.fq>

**Post-alignment**

7. Convert reads to BAM and sort reads by position with [samtools](https://github.com/samtools/samtools).

This pipeline is just a template; it can be adjusted easily in code. Each step
is an ordered dictionary (to preserve order), where keys indicate which path to
take. So far, only the pre-processing steps are configurable.

The pipeline implements two steps (which are hidden from the user):

First, the *dispatch* step creates a text file of all the sample's run
information (called the `<job>_samples.txt` file). Each of these is a JSON
entry in a text file, which is passed to a runner. Each of these JSON entries
contains all data needed by the processing command: reference, parameters,
sample information, etc.

Then, the dispatch command writes a Slurm batch script in your directory for
this job, and runs it. This batch script then calls the *runners* (using
Slurm's array job infrastructure) which is a subcommand of this pipeline that
takes a single line from the `<job>_samples.txt` file and run it.

## Running the Pre-Processing and Alignment Pipeline

All processing is done on a single pair of reads. In order to run, two
configuration files are needed (note you can name these whatever you want):

1. A `setup.json` JSON configuration file, which includes paths to all programs.

2. A `samples.txt` tab-delimited file which provides a mapping between both
	 read pairs, sample IDs (`@RG:SM` in the SAM spec), and read group IDs
   (`@RG:ID` in the SAM spec). The column format is:

        sample_id, read_id, reads1, reads2

You'll need to have have Python installed with
[logbook](https://github.com/mitsuhiko/logbook). Logbook is used to log
info/errors. Note that this uses Slurm's logging mechanism; task ID gets it's
own log. We might add ZeroMQ-based messaging for log consolidation. On Slurm,
tell load Python with:

    $ module load python

Then start you run with something like:

    $ python sample_dispatch.py dispatch --ref ~/data/refgen3-hapmap/maize3.fa.gz \
      --scythe --pre-seqqs --post-seqqs --trimfq \
      --samples samples.txt --adapter illumina_adapters.fa --log log --stats data/stats/ \
      --mem 1024 --threads 10 --job teoparent --setup setup.json  --bam-dir data/aln/ \
      -P bigmem

This will create `<job>_samples.txt` and `<job>_batch.sh`. The first file is
the job sample file, the second file is the Slurm batch script. You do *not*
need to start the batch script; this is just for reproducibility.
`sample_dispatch.py` will automatically start the batch script and notify you
if it has been submitted. You can prevent this by specify `--dry-run`.

## Todo

This software is in alpha. Some serious limitations are:

- Logging is implemented with Slurm task IDs; there is no consolidated logging
  across array jobs.

- Fragile exception raising during validation; the validation is a bit crufty.

- Commands are not modular.

- The `setup.json` file doesn't have a method to handle modules. We can't count
  on trying `module load` on each program, as modules can load many programs
  (with different paths).



