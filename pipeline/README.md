

This document provides a general workflow and overview of the tools we
have used to analyse rRNA processing and modification in Archaea using
Nanopore RNA-seq data, including:  
- Basecalling and demultiplexing of raw FAST5 reads using `guppy`  
- Trimming of reads using `pychopper` & `cutadapt`  
- Mapping of reads to the genome using `minimap2`  
- Read coverage analysis using `samtools`  
- Modified base detection using `tombo`, `Eligos2` & `pysamstats`

------------------------------------------------------------------------

#### Table of Contents

-   <a href="#library-preparation" id="toc-library-preparation">Library
    preparation</a>
    -   <a href="#direct-rna" id="toc-direct-rna">Direct RNA</a>
    -   <a href="#direct-cdna" id="toc-direct-cdna">Direct cDNA</a>
-   <a href="#sequencing" id="toc-sequencing">Sequencing</a>
-   <a href="#data-analysis" id="toc-data-analysis">Data analysis</a>
    -   <a href="#data-management" id="toc-data-management">Data management</a>
    -   <a
        href="#basecalling-demultiplexing-and-trimming-of-direct-rna-libraries"
        id="toc-basecalling-demultiplexing-and-trimming-of-direct-rna-libraries">Basecalling,
        demultiplexing and trimming of direct RNA libraries</a>
        -   <a href="#demulitplexing-using-poreplex"
            id="toc-demulitplexing-using-poreplex">Demulitplexing using
            <span><code>poreplex</code></span></a>
        -   <a href="#basecalling-using-guppy"
            id="toc-basecalling-using-guppy">Basecalling using
            <code>guppy</code></a>
        -   <a href="#polya-trimming-using-cutadapt"
            id="toc-polya-trimming-using-cutadapt">Poly(A)-trimming using
            <span>cutadapt</span></a>
    -   <a
        href="#basecalling-demultiplexing-and-trimming-of-direct-cdna-libraries"
        id="toc-basecalling-demultiplexing-and-trimming-of-direct-cdna-libraries">Basecalling,
        demultiplexing and trimming of direct cDNA libraries</a>
        -   <a href="#basecalling-using-guppy-1"
            id="toc-basecalling-using-guppy-1">Basecalling using
            <code>guppy</code></a>
        -   <a href="#demultiplexing-of-basecalled-reads-using-guppy_barcoder"
            id="toc-demultiplexing-of-basecalled-reads-using-guppy_barcoder">Demultiplexing
            of basecalled reads using <code>guppy_barcoder</code></a>
        -   <a
            href="#read-orientation-and-detection-of-full-length-sequenced-reads-using-pychopper"
            id="toc-read-orientation-and-detection-of-full-length-sequenced-reads-using-pychopper">Read
            orientation and detection of full-length sequenced reads using
            <span><code>pychopper</code></span></a>
    -   <a href="#read-alignment-using-minimap2"
        id="toc-read-alignment-using-minimap2">Read alignment using
        <span><code>minimap2</code></span></a>
    -   <a href="#calculation-of-coverage-files"
        id="toc-calculation-of-coverage-files">Calculation of coverage files</a>
    -   <a
        href="#detection-of-rrna-processing-sites-and-classification-of-rrna-intermediates"
        id="toc-detection-of-rrna-processing-sites-and-classification-of-rrna-intermediates">Detection
        of rRNA processing sites and classification of rRNA intermediates</a>
    -   <a href="#circular-rna-detection"
        id="toc-circular-rna-detection">Circular RNA detection</a>
    -   <a href="#modified-base-detection"
        id="toc-modified-base-detection">Modified base detection</a>

------------------------------------------------------------------------

## Library preparation

### Direct RNA

Libraries for Nanopore direct RNA sequencing (*DRS*) were prepared from
poly(A)-tailed RNAs according to the SQK-RNA001 Kit protocol (Oxford
Nanopore, Version: DRS_9026_v1_revP_15Dec2016) with minor modifications.
Detailed protocols can be found in the [Oxford
Nanopore](https://nanoporetech.com) community.

### Direct cDNA

Libraries for direct cDNA sequencing (*cDNA*) were prepared following
the instructions in the direct cDNA sequencing with native barcoding
protocol (SQK-DCS109 with EXP-NBD104) from Oxford Nanopore Technologies
with minor modifications. Briefly, the VN Primer was replaced with a
custom 3’ cDNA RT primer (5′-
ACTTGCCTGTCGCTCTATCTTCATTGATGGTGCCTACAG-3′, 2 µM).

## Sequencing

*DRS* and *cDNA* libraries were sequenced on a MinION Mk1B or Mk1C using
R9.4 flow cells and the recommended scripts in MinKNOW to generate FAST5
files.

> Note: Live-basecalling in fast mode was enabled to monitor
> translocation speed and quality during a run.

## Data analysis

### Data management

To simplify readability of the following scripts and the Rscripts, the
folder management is shown in the following (only exemplarily for the
direct cDNA datasets, direct RNA paths are commented when necessary):

``` bash
rRNA_maturation/
└── MinKNOW_output
└── rebasecallling/
    ├── sequencing_summary.txt
    ├── pass
    └── fail
└── demultiplexing/
    ├── barcoding_summary.txt
    ├── pass
    └── fail
└── analysis/
    ├── fastq_pass
    ├── fastq_pass_fl
    ├── summary
    ├── genome
    ├── sample_info
    ├── pychopper_edlib
    ├── pychopper_edlib_rescue
    ├── mapped_fl_only
    ├── coverage
    └── coverage_fl_only
```

### Basecalling, demultiplexing and trimming of direct RNA libraries

#### Demulitplexing using [`poreplex`](https://github.com/hyeshik/poreplex)

Before starting a sequencing run, different running options can be
selected in MinKNOW. In case reads are stored in multi-FAST5-containing
files (e.g. 4000 reads per file), files can be converted to single-read
FAST5 files using the
[ont_fast5_api](https://github.com/nanoporetech/ont_fast5_api), as some
workflows
(e.g. [`nanopolish`](https://nanopolish.readthedocs.io/en/latest/) and
[`tombo`](https://nanoporetech.github.io/tombo/)) rely on single-FAST5
files for further analysis.  
After a run, reads are stored in two folders (*fast5_failed*,
*fast5_passed*). To prevent actual good reads from beeing discarded we
**included all reads from both folders** in the following steps of the
analysis.  
First, we converted multi-FAST5-files with the `multi_to_single_fast5`
command from the
[ont_fast5_api](https://github.com/nanoporetech/ont_fast5_api):

``` bash
# set input_path, save_path, search in all folders for fast5 files, set number of threads
multi_to_single_fast5 \
    --input_path \   # path folder containing multi_read_fast5 files  
    --save_path \    # path to folder where single_read fast5 files will be output  
    --recursive \    # recursively search sub-directories  
    --threads \      # number of CPU threads to use  
```

The output will be single-read FAST5 files in the *save_path* folder
with one subfolder per multi-read input file.

Multiplexed libraries (*how to* is described here:
<https://github.com/hyeshik/poreplex>) can be demultiplexed using
[`poreplex`](https://github.com/hyeshik/poreplex). Following this
approach four direct RNA sequencing libraries can be barcoded, pooled
and sequenced together. `Poreplex` can demultiplex the libraries into
separate folders with:

``` bash
# trim adapters, basecall using albacore, de-multiplex, create symbolic fast5 link, increase number of working processes, sort reads to folders according to barcodes
poreplex \
    -i \                 # path/to/fast5   
    -o \                 # path/to/output
    --trim-adapter \     # trim 3′ adapter sequences from FASTQ outputs   
    --barcoding \        # sort barcoded reads into separate outputs
    --basecall  \        # call the ONT albacore for basecalling on-the-fly  
    --symlink-fast5 \    # create symbolic links to FAST5 files in output directories even when hard linking is possible
    --parallel \         # number of worker processes
```

Reads can be basecalled automatically during demultiplexing using
`albacore`, the outdated basecaller of ONT. As `albacore` is meanwhile
outdated for a long time, we chose to only sort the reads based on
`poreplex` and used `guppy` (Version 6.1.3+cc1d765d3) for basecalling in
high-accuracy mode. In this case the `poreplex` can be shortened to:
`poreplex -i <input> -o <output> --barcoding --parallel`.

#### Basecalling using `guppy`

Demultiplexed raw FAST5 files can be basecalled (*translated* in FASTQ
data) using `guppy`, the ONT-developed basecaller (available in the ONT
Community). We used Version 6.1.3+cc1d765d3 for basecalling of all of
our reads:

``` bash
guppy_basecaller \
--input_path ${input} \          # input path    
--save_path ${output} \          # output path      
-c rna_r9.4.1_70bps_hac.cfg  \   # high-accuracy DRS config
--compress_fastq \
--fast5_out \
--recursive \
--progress_stats_frequency 60 \
-x 'auto'
```

#### Poly(A)-trimming using [cutadapt](https://cutadapt.readthedocs.io/en/stable/)

After basecalling and demultiplexing, artificially added polyA sequences
were removed from the 3’ ends using cutadapt v4.1 (-a A{10}, -e 3, -j 0)
X.

``` bash
for file in ${ext_dir}/rebasecalling/*/*.fq.gz
do 
  filename_extended=${file##*/}
  filename=${filename_extended%%.*}
  foldername=$(echo $file | cut -d"/" -f 6)
  foldername2=$(echo $filename | rev | cut -d"_" -f 1 | rev)

  out_py=${ext_dir}/analysis/cutadapt_fastq/${foldername}
  mkdir -p ${out_py}

  cutadapt \
    -a "A{10}" \
    -e 3 \
    -j 0 \
    -o ${out_py}/${foldername}_${foldername2}_cutadapt.fastq \
    ${file}
done

# Merge all (WT haloferax and ∆KsgA haloferax)
cat ${ext_dir}/analysis/cutadapt_fastq/hvo_notex/*.fastq > ${ext_dir}/analysis/cutadapt_fastq/hvo_notex_cut.fastq
cat ${ext_dir}/analysis/cutadapt_fastq/hvo_notex_dksga/*.fastq > ${ext_dir}/analysis/cutadapt_fastq/hvo_notex_dksga_cut.fastq
```

### Basecalling, demultiplexing and trimming of direct cDNA libraries

#### Basecalling using `guppy`

After sequencing (and despite live-basecalling) all datasets in the
raw_FAST5 📁 were re-basecalled using `guppy` (v. 6.3.2+bb5453e) in
high-accuracy mode with a q-score cutoff of 7. Next, basecalled files
were demultiplexed using the `guppy_barcoder` command from the `guppy`
suite (available in the [ONT community](https://nanoporetech.com) using
default parameters). Demultiplexed files in FASTQ format were written to
the fastq_pass 📁.

``` bash
# files
input=rRNA_maturation/MinKNOW_output
output_cDNA=rRNA_maturation/rebasecalling

# Basecalling of cDNA files 
guppy_basecaller \
--input_path ${input} \          # input path
--save_path ${output_cDNA} \     # output path
-c dna_r9.4.1_450bps_hac.cfg \   # config file: high accuracy cDNA 
--compress_fastq \               # compress output
--fast5_out \                    # output FAST5
--min_qscore 7 \                 # qscore cutoff
--recursive \                    # look for FAST5 recursively in path
--progress_stats_frequency 60 \  # output progress every minute
-x auto                          # automatic detection of GPUs
```

> DRS & cDNA runs require different options.  
> Config file selection based on selected accuracy, flowcell version,
> library preparation kit are listed with
> `guppy_basecaller --print_workflows`

Using the selected options `guppy` produces fast5_pass, fast5_fail,
fastq_pass, fastq_fail, summary and report files that are written to the
rebasecallling 📁. Multiple FASTQs can be merged using
`cat rRNA_maturation/rebasecallling/*.fastq > rRNA_maturation/rebasecallling/run_id.fastq`.

Sequencing summary files are also written to the rebasecallling 📂 and
are used during the quality control of the runs and reads. For better
viewing they can be moved to the analysis/summary 📁 using
`mv rRNA_maturation/rebasecalling/sequencing_summary.txt rRNA_maturation/analysis/summary/sequencing_summary.txt`

#### Demultiplexing of basecalled reads using `guppy_barcoder`

Next, multiplexed cDNA libraries are demultiplexed in a separate step
using `guppy_barcoder` (v. 6.3.2+bb5453e) in default mode (barcode kit
EXP-NBD104).

``` bash
# files
input=rRNA_maturation/rebasecallling
output=rRNA_maturation/demultiplexing

# Demultiplexing of cDNA files   
## fail folder 
guppy_barcoder \
--input_path $input/fail \
--save_path $output/fail \
--recursive \
--progress_stats_frequency 60 \
-x 'auto' \
--compress_fastq \
--barcode_kits 'EXP-NBD104' 

## pass folder 
guppy_barcoder \
--input_path $input/pass \
--save_path $output/pass \
--recursive \
--progress_stats_frequency 60 \
-x 'auto' \
--compress_fastq \
--barcode_kits 'EXP-NBD104' 
```

Multiple FASTQs are written to the demultiplexing 📂 and can be merged
with
e.g. `cat rRNA_maturation/demultiplexing/pass/barcode01/*.fastq > rRNA_maturation/analysis/fastq_pass/barcode01.fastq`.
Barcode summary files are written to the rRNA_maturation/demultiplexing
📁 and can be moved to the rRNA_maturation/analysis/summary 📂.

#### Read orientation and detection of full-length sequenced reads using [`pychopper`](https://github.com/nanoporetech/pychopper)

Full-length cDNA reads containing SSP (TTTCTGTTGGTGCTGATATTGCTGGG) and
custom VNP primers (ACTTGCCTGTCGCTCTATCTTCATTGATGGTGCCTACAG) in the
correct orientation were identified using `pychopper` (v.2.5.0) with
standard parameters using the edlib backend and autotuned cutoff
parameters estimated from subsampled data.  
\> Note that a custom primer fasta file was created as described in
<https://github.com/epi2me-labs/pychopper>

Additionally, protocol-specific read rescue was performed using the
DCS109-specific command in pychopper. All full-length detected reads
were merged and used for subsequent steps. The output is saved in
pychopper_edlib 📁.

``` bash
# files
ext_dir=rRNA_maturation

# perform pychopper for all cDNA files
for file in ${ext_dir}/analysis/fastq_pass/*/*.fq
do 
  filename_extended=${file##*/}
  filename=${filename_extended%%.*}

  echo ${filename} pychopper edlib started ...

  out_py=${ext_dir}/analysis/pychopper_edlib/${filename}
  mkdir -p ${out_py}
  
  pychopper \
  -m edlib \
  -b pychopper_custom.fasta \
  -r ${out_py}/${filename}_report.pdf \
  -t 8 \
  -u ${out_py}/${filename}_unclassified.fastq \
  -w ${out_py}/${filename}_rescued.fastq \
  -S ${out_py}/${filename}_stats.txt \
  ${file} \
  ${out_py}/${filename}_full_length.fastq

  # After a first round, a second round of `pychopper` was applied to the unclassified direct cDNA reads with DCS-specific read rescue enabled. 
  echo ${filename} DCS109 specific rescue started ...
  out_py2=${ext_dir}/analysis/pychopper_edlib_rescue/${filename}
  mkdir -p ${out_py2}
  
  pychopper \
  -x DCS109 \
  -b pychopper_custom.fasta \
  -r ${out_py2}/${filename}_report.pdf \
  -t 8 \
  -m edlib \
  -u ${out_py2}/${filename}_unclassified.fastq \
  -w ${out_py2}/${filename}_rescued.fastq \
  -S ${out_py2}/${filename}_stats.txt \
  ${out_py}/${filename}_unclassified.fastq \
  ${out_py2}/${filename}_full_length.fastq
  
  echo ${filename} files are merged ...
  out_py3=${ext_dir}/analysis/fastq_pass_fl/${filename}
  mkdir -p ${out_py3}
  
  # finally, all fl files are merged
  cat ${out_py}/${filename}_full_length.fastq ${out_py2}/${filename}_full_length.fastq > ${out_py3}/${filename}_full_length.fastq

done
```

### Read alignment using [`minimap2`](https://github.com/lh3/minimap2)

Direct RNA and direct cDNA reads were mapped using minimap2 (v.
2.24-r1122) with standard parameters suggested for the alignment of
noisy direct RNA reads (-ax splice -uf -k14) and Nanopore genomic reads
(-max map-ont), respectively. In each case, –MD tag was used to include
the MD tag for calculating mapping identities. Alignment files were
converted to BAM files, sorted and indexed using samtools (v.1.15.1).

Files were mapped to the representative reference genomes (downloaded to
genome 📁) from  
- [*Haloferax
volcanii*](https://www.ncbi.nlm.nih.gov/genome/?term=haloferax+volcanii)  
- [*Sulfolobus
acidocaldarius*](https://www.ncbi.nlm.nih.gov/genome/?term=Sulfolobus+acidocaldarius) -
[*Pyrococcus
furiosus*](https://www.ncbi.nlm.nih.gov/genome/?term=Pyrococcus+furiosus)

and to circularly permuted 16S and 23S rRNA sequences (which will be
described in on of the next sections)

To facilitate mapping, a sample information file was used to
automatically switch between reference genomes based on the barcode id.

``` r
# R!  

# load libraries ----
library(vroom)
library(data.table)

# make table ----
sample_info <- data.table::data.table(bc  = sprintf("barcode%02d", c(1,2,4,5,7,8)),
                                      org = c(rep("Haloferax volcanii H26",2),
                                              rep("Sulfolobus acidocaldarius",2),
                                              rep("Pyrococcus furiosus",1),
                                              rep("Haloferax volcanii H26",1)),
                                      org_s = c(rep("hvo",2),
                                                rep("sac",2),
                                                rep("pfu",1),
                                                rep("hvo",1)),
                                      org_z = c(rep("hvo.transcripts",2),
                                                rep("sac.transcripts",2),
                                                rep("pfu.transcripts",1),
                                                rep("hvo.transcripts",1)),
                                      cond = c("OD02_rep1",
                                               "OD02_rep2",
                                               "OD02_rep1",
                                               "OD02_rep2",
                                               "95_rep1",
                                               "IVT"))

ext_dir <- "rRNA_maturation"
vroom_write(sample_info, paste0(ext_dir, "/analysis/sample_info/sample_info.txt"), delim = ",", col_names = F)
```

Finally, reads were aligned to the reference genomes using:

``` bash
ext_dir="rRNA_maturation"

# direct cDNA 
for file in ${ext_dir}/analysis/fastq_pass_fl/*/*.fastq
do 
  filename_extended=${file##*/}
  filename=${filename_extended%%.*}
  foldername=$(echo $filename_extended | cut -d"_" -f 1 | cut -d"." -f 1)
  
  echo ${foldername} mapping started ...
  
  out_py=${ext_dir}/analysis/mapped_fl_only/${foldername}
  mkdir -p ${out_py}
  
  map=$(grep "$foldername" ${ext_dir}/analysis/sample_info/sample_info.txt | cut -d',' -f 3)

  if [[ $map != "NA" ]]; 
    then
      echo $map
      minimap2 -ax map-ont --MD ${ext_dir}/analysis/genome/${map}.fasta $file > ${out_py}/${foldername}".sam"
      samtools flagstat ${out_py}/${foldername}".sam" > ${out_py}/${foldername}"_stats.txt"
      samtools sort ${out_py}/${foldername}".sam" -o ${out_py}/${foldername}".sorted.bam"
      samtools index ${out_py}/${foldername}".sorted.bam"
      rm ${out_py}/${foldername}".sam"
  fi
done

# direct RNA (RNA files saved to similar folder structure in different project)
## use polyA-trimmed files 
for file in ${ext_dir}/analysis/cutadapt_fastq/*_cut.fastq
do 
  filename_extended=${file##*/}
  filename=${filename_extended%%.*}
  fix=$(echo $filename | cut -d"c" -f 1)
  foldername=${fix%?}
  
  echo ${foldername} mapping started ...

  out_py=${ext_dir}/analysis/mapped_fl_only/${foldername}
  mkdir -p ${out_py}
  
  minimap2 -ax splice -uf -k14 --MD ${ext_dir}/analysis/genome/hvo.fasta $file > ${out_py}/${foldername}".sam"
  samtools flagstat ${out_py}/${foldername}".sam" > ${out_py}/${foldername}"_stats.txt"
  samtools view -bS ${out_py}/${foldername}".sam" -o ${out_py}/${foldername}".bam"
  samtools sort ${out_py}/${foldername}".bam" -o ${out_py}/${foldername}".sorted.bam"
  samtools index ${out_py}/${foldername}".sorted.bam"
  rm ${out_py}/${foldername}".sam"
  rm ${out_py}/${foldername}".bam"
done
```

### Calculation of coverage files

Coverage profiles (compare Figure 1C, Supplementary Figure 4) were
calculated using
[`samtools depth`](http://www.htslib.org/doc/samtools-depth.html) with
the options:  
- a (*to output all positions*)  
- J (*to include reads with deletions in depth calculations*)

``` bash
ext_dir="rRNA_maturation"

for file in ${ext_dir}/analysis/fastq_pass_fl/*/*.fastq
do 
  filename_extended=${file##*/}
  filename=${filename_extended%%.*}
  foldername=$(echo $filename_extended | cut -d"." -f 1)
  
  echo ${foldername} coverage calc started ...
  
  # > all passed
  out_py=${ext_dir}/analysis/coverage/${foldername}
  mkdir -p ${out_py}
  
  samtools depth -a -J ${ext_dir}/analysis/mapped/${foldername}/${foldername}".sorted.bam" -o ${out_py}/${foldername}_cov.tsv
  
  # > full-length
  out_py2=${ext_dir}/analysis/coverage_fl_only/${foldername}
  mkdir -p ${out_py2}
  
  samtools depth -a -J ${ext_dir}/analysis/mapped_fl_only/${foldername}/${foldername}".sorted.bam" -o ${out_py2}/${foldername}_cov.tsv
done
```

### Detection of rRNA processing sites and classification of rRNA intermediates

### Circular RNA detection

### Modified base detection
