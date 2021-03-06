
# MNase-seq analysis pipeline

## description

An analysis pipeline for paired-end MNase-seq data with the following major steps:

- quality trimming with [cutadapt](http://cutadapt.readthedocs.io/en/stable/guide.html)
- alignment with [bowtie1](http://bowtie-bio.sourceforge.net/index.shtml)
- selection of correctly paired reads
- summaries of quality statistics from [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 
- summaries of library processing statistics
- summaries of library fragment sizes
- generation of nucleosome dyad (i.e. fragment midpoint) and nucleosome protection coverage tracks
- library size and spike-in normalization of coverage
- genome-wide scatterplots and correlations
- quantification and visualization of nucleosome occupancy, fuzziness, and position shifts with [DANPOS2](https://sites.google.com/site/danposdoc/)
- differential occupancy analysis over transcripts with [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- data visualization (heatmaps and metagenes, with the option to separate data into clusters of similar signal)

## requirements

### required software

- Unix-like operating system (tested on CentOS 7.2.1511)
- Git
- [conda](https://conda.io/docs/user-guide/install/index.html)
- [build-annotations pipeline](https://github.com/winston-lab/build-annotations)

### required files

- Paired-end FASTQ files of MNase-seq libraries prepared as described in [our publication](https://doi.org/10.1016/j.molcel.2018.09.005). FASTQ files should be demultiplexed, with a separate file for read 1 and read2 and all 5' inline barcodes trimmed. A separate pipeline for demultiplexing paired-end FASTQ files with 5' inline barcodes can be found [here](https://github.com/winston-lab/demultiplex-paired-end). This pipeline has only been tested using Illumina sequencing data. 

- FASTA files:
    - the 'experimental' genome
    - the spikein genome, if any samples have spikeins

- [BED6](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format annotation files:
    - transcript annotation
    - optional: other annotations for data visualization (i.e. heatmaps and metagenes)

## instructions

**0**. If you need to demultiplex and trim your FASTQ files, use the separate ['demultiplex-paired-end' pipeline](https://github.com/winston-lab/demultiplex-paired-end) to do so.

**1**. If you haven't already done so, clone the separate ['build-annotations' pipeline](https://github.com/winston-lab/build-annotations), make a copy of the `config_template.yaml` file called `config.yaml`, and edit `config.yaml` as needed so that it points to the experimental genome FASTA file and a transcript annotation. Motif database building and paths to other annotations used to create region annotation files can be disabled if you are not using another pipelines that uses these attributes.

```bash

# clone the repository
git clone https://github.com/winston-lab/build-annotations.git

# move into the build-annotations pipeline directory
cd build-annotations

# make a copy of the configuration template file
cp config_template.yaml config.yaml

# edit the configuration file
vim config.yaml         # or use your favorite editor
```

**2**. Clone this repository.

```bash
git clone https://github.com/winston-lab/mnase-seq.git
```

**3**. Create and activate the `snakemake_default` virtual environment for the pipeline using conda. The virtual environment creation can take a while. If you've already created the `snakemake_default` environment from another one of my pipelines, this is the same environment, so you can skip creating the environment and just activate it.

```bash
# navigate into the pipeline directory
cd mnase-seq

# create the snakemake_default environment
conda env create -v -f envs/snakemake_default.yaml

# activate the environment
source activate snakemake_default

# to deactivate the environment
# source deactivate
```

**4**. Make a copy of the configuration file template `config_template.yaml` called `config.yaml`, and edit `config.yaml` to suit your needs.

```bash
# make a copy of the configuration template file
cp config_template.yaml config.yaml

# edit the configuration file
vim config.yaml    # or use your favorite editor
```

**5**. With the `snakemake_default` environment activated, do a dry run of the pipeline to see what files will be created.

```bash
snakemake -p --use-conda --dryrun
```

**6**. If running the pipeline on a local machine, you can run the pipeline using the above command, omitting the `--dryrun` flag. You can also use N cores by specifying the `--cores N` flag. The first time the pipeline is run, conda will create separate virtual environments for some of the jobs to operate in. Running the pipeline on a local machine can take a long time, especially for many samples, so it's recommended to use an HPC cluster if possible. On the HMS O2 cluster, which uses the SLURM job scheduler, entering `sbatch slurm_submit.sh` will submit the pipeline as a single job which spawns individual subjobs as necessary. This can be adapted to other job schedulers and clusters by adapting `slurm_submit.sh`, which submits the pipeline to the cluster, `slurm_status.sh`, which handles detection of job status on the cluster, and `cluster.yaml`, which specifies the resource requests for each type of job.

