

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

- <a href="#library-preparation" id="toc-library-preparation">Library
  preparation</a>
  - <a href="#direct-rna" id="toc-direct-rna">Direct RNA</a>
  - <a href="#direct-cdna" id="toc-direct-cdna">Direct cDNA</a>
- <a href="#sequencing" id="toc-sequencing">Sequencing</a>
- <a href="#data-analysis" id="toc-data-analysis">Data analysis</a>
  - <a href="#data-management" id="toc-data-management">Data management</a>
  - <a
    href="#basecalling-demultiplexing-and-trimming-of-direct-rna-libraries"
    id="toc-basecalling-demultiplexing-and-trimming-of-direct-rna-libraries">Basecalling,
    demultiplexing and trimming of direct RNA libraries</a>
    - <a href="#demulitplexing-using-poreplex"
      id="toc-demulitplexing-using-poreplex">Demulitplexing using
      <code>poreplex</code></a>
    - <a href="#basecalling-using-guppy"
      id="toc-basecalling-using-guppy">Basecalling using
      <code>guppy</code></a>
    - <a href="#polya-trimming-using-cutadapt"
      id="toc-polya-trimming-using-cutadapt">Poly(A)-trimming using
      cutadapt</a>
  - <a
    href="#basecalling-demultiplexing-and-trimming-of-direct-cdna-libraries"
    id="toc-basecalling-demultiplexing-and-trimming-of-direct-cdna-libraries">Basecalling,
    demultiplexing and trimming of direct cDNA libraries</a>
    - <a href="#basecalling-using-guppy-1"
      id="toc-basecalling-using-guppy-1">Basecalling using
      <code>guppy</code></a>
    - <a href="#demultiplexing-of-basecalled-reads-using-guppy_barcoder"
      id="toc-demultiplexing-of-basecalled-reads-using-guppy_barcoder">Demultiplexing
      of basecalled reads using <code>guppy_barcoder</code></a>
    - <a
      href="#read-orientation-and-detection-of-full-length-sequenced-reads-using-pychopper"
      id="toc-read-orientation-and-detection-of-full-length-sequenced-reads-using-pychopper">Read
      orientation and detection of full-length sequenced reads using
      <code>pychopper</code></a>
  - <a href="#read-alignment-using-minimap2"
    id="toc-read-alignment-using-minimap2">Read alignment using
    <code>minimap2</code></a>
  - <a href="#analysis-of-coverage-files"
    id="toc-analysis-of-coverage-files">Analysis of coverage files</a>
    - <a href="#calculation" id="toc-calculation">Calculation</a>
    - <a href="#plotting" id="toc-plotting">Plotting</a>
  - <a href="#quality-control" id="toc-quality-control">Quality control</a>
  - <a
    href="#detection-of-rrna-processing-sites-and-classification-of-rrna-intermediates"
    id="toc-detection-of-rrna-processing-sites-and-classification-of-rrna-intermediates">Detection
    of rRNA processing sites and classification of rRNA intermediates</a>
  - <a href="#circular-rna-detection"
    id="toc-circular-rna-detection">Circular RNA detection</a>
  - <a href="#modified-base-detection"
    id="toc-modified-base-detection">Modified base detection</a>
    - <a href="#stage-processing" id="toc-stage-processing">Stage
      processing</a>
    - <a href="#stage-sorting" id="toc-stage-sorting">Stage sorting</a>
    - <a href="#esb-calculation" id="toc-esb-calculation">ESB calculation</a>
    - <a href="#eligos2" id="toc-eligos2">Eligos2</a>
    - <a href="#modified-base-detection-based-on-signal-data"
      id="toc-modified-base-detection-based-on-signal-data">Modified base
      detection based on signal data</a>

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
the fastq_pass 📂.

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

Sequencing summary files are also written to the rebasecallling 📁 and
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
genome 📂) from  
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

### Analysis of coverage files

#### Calculation

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

#### Plotting

``` r
### R ###

# libraries ---- 
library(here)
library(vroom)
library(tidyverse)
library(ape)

# functions ----
read_in_gff_rrna <- function(input_file){
  read.gff(input_file) %>%
    dplyr::filter(type %in% "rRNA") %>%
    as_tibble() %>%
    dplyr::mutate(start_feature = start, end_feature = end,strand_feature = strand) %>%
    dplyr::mutate(Parent = str_split_fixed(str_split_fixed(attributes, ";Parent=",2)[,2],";Dbxref",2)[,1],
                  ecogene = str_split_fixed(str_split_fixed(attributes, ",GeneID", 2)[,1], "EcoGene:",2)[,2],
                  short_gene = str_split_fixed(str_split_fixed(attributes, ";locus_tag", 2)[,1], "gene=",2)[,2],
                  id_name = str_split_fixed(str_split_fixed(attributes, ";Parent=", 2)[,1], "ID=", 2)[,2],
                  locus_name = str_split_fixed(str_split_fixed(attributes, ";product=", 2)[,2], " ", 2)[,1], 
                  width = abs(start_feature - end_feature)) %>%
    dplyr::select(seqid, id_name, locus_name, start_feature, end_feature, strand_feature, Parent, type, width)
}

read_in_cov <- function(input,myBC, myMet){
  vroom(input, col_names = c("chr", "pos", "cov")) %>%
    dplyr::mutate(barcode = myBC,
                  method = myMet,
                  max_cov = max(cov,na.rm = T),
                  cov_norm = cov/max_cov)
}

get_org_filter <- function(input){
  input %>%
    left_join(org_t, by = "chr") %>%
    group_by(org) %>%
    dplyr::filter(pos %in% startpos:endpos)
}

# data ----
## genome data ====
### fasta ####
hvo_fasta <- readDNAStringSet(here("data/genome/hvo.fasta"))
names(hvo_fasta) <- str_split_fixed(names(hvo_fasta), " ",2)[,1]
sac_fasta <- readDNAStringSet(here("data/genome/sac.fasta"))
names(sac_fasta) <- str_split_fixed(names(sac_fasta), " ",2)[,1]
pfu_fasta <- readDNAStringSet(here("data/genome/pfu.fasta"))
names(pfu_fasta) <- str_split_fixed(names(pfu_fasta), " ",2)[,1]

### gff ####
hvo_gff_rrna <- read_in_gff_rrna(here("data/genome/hvo.gff"))
sac_gff_rrna <- read_in_gff_rrna(here("data/genome/sac.gff"))
pfu_gff_rrna <- read_in_gff_rrna(here("data/genome/pfu.gff"))

### organism rRNA coordinate table ####
org_t <- data.table(chr = c(names(hvo_fasta),names(sac_fasta), names(pfu_fasta)),
                    org = c(rep("hvo", 5), "sac", "pfu"),
                    startpos = c(rep(hvo_gff_rrna$start_feature[1]-500,5),
                                 rep(sac_gff_rrna$start_feature[1]-500,1),
                                 rep(pfu_gff_rrna$start_feature[1]-500,1)),
                    endpos = c(rep(hvo_gff_rrna$end_feature[3]+500,5),
                               rep(sac_gff_rrna$end_feature[2]+500,1),
                               rep(pfu_gff_rrna$end_feature[2]+500,1)))
## all mapped ====
cov_files <- list.files(path = paste0(dir,"coverage_all"), pattern = "_cov.tsv", recursive = T, full.names = T)
cov_bc    <- str_extract(cov_files, "barcode\\d+")
mapped_all <- pmap_dfr(list(cov_files[c(1:9)], cov_bc[c(1:9)], "all"), read_in_cov)
mapped_all_f <- get_org_filter(mapped_all)

## full-length ====
### all full-length ####
cov_fl_files <- list.files(path = paste0(dir,"coverage_fl_only"), pattern = "_cov.tsv", recursive = T, full.names = T)
fl_bc       <- str_extract(cov_fl_files, "barcode\\d+")
mapped_fl <- pmap_dfr(list(cov_fl_files[c(1:9)], fl_bc[c(1:9)], "fl"), read_in_cov)
mapped_fl_f <- get_org_filter(mapped_fl)

## combine ====
mapped_cov <- bind_rows(mapped_all_f, mapped_fl_f)

# plotting ----
## e.g. for barcode 01 (Figure 1C) ====
ggplot(data = mapped_cov %>% 
         dplyr::filter(method %in% c("all"),
                       barcode %in% "barcode01"),
       aes(x = pos, y = cov, fill = method)) +
  geom_area(position = "identity", alpha = 0.8) +
  geom_line(color = "black") +
  theme_pubclean() +
  theme(panel.grid.minor = element_blank()) +
  ylab("Nr of reads at position") +
  xlab("") +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = c("#99A8D8", "#001959")) 
```

### Quality control

### Detection of rRNA processing sites and classification of rRNA intermediates

We used custom R scripts to evaluate accuracy of known terminal rRNA
positions, including mature positions and bulge-helix-bulge-cleavage
sites.

For rRNA stage classification, coordinates of 5’ and 3’ ends, clipping
information extracted from the CIGAR string, strand information and
sequence identity were considered. Briefly, only reads with a mapping
identity ≥ 90 % and with soft-clippings ≤ 20 nucleotides were used for
classification based on linear templates (reads mapping to the
representative genomes). Finally, stage-sorted reads were plotted in a
genome browser-like view and evaluated based on the terminal positions.
Enriched and previously undescribed combinations of connected 5’ and 3’
ends were considered in the classification.

``` r
### R ###

# libraries ----
source(here("Rscripts/load_libraries.R"))


# functions ----
read_in_bam_large <- function(file){
  
  init <- readGAlignments(file, use.names = T, param = ScanBamParam(tag=c("NM"), what="mapq"))
  init_t <- GenomicAlignments::as.data.frame(init) %>%
    dplyr::mutate(minion_read_name = names(init),
                  mapped_gene = seqnames) 
  
  left  <- paste(str_split_fixed(string = init_t$cigar, pattern = "M", n = 2)[,1],"M", sep = "")
  right <- paste(str_split_fixed(string = init_t$cigar, pattern = "M", n = 2)[,2],"1M", sep = "")
  
  #................................calculate cigar tables / SOFT AND HARD CLIPPING!!!
  init_t$soft_l <- as_tibble(cigarOpTable(left))$S
  init_t$hard_l <- as_tibble(cigarOpTable(left))$H
  init_t$soft_r <- as_tibble(cigarOpTable(right))$S
  init_t$hard_r <- as_tibble(cigarOpTable(right))$H
  init_t$matches <- as_tibble(cigarOpTable(init_t$cigar))$M
  init_t$insertions <- as_tibble(cigarOpTable(init_t$cigar))$I
  init_t$deletions <- as_tibble(cigarOpTable(init_t$cigar))$D
  init_t$skip <- as_tibble(cigarOpTable(init_t$cigar))$N
  init_t$cigar_S <- as_tibble(cigarOpTable(init_t$cigar))$S
  
  # calculate number of aligned reads based on CIGAR operations (M,I)
  init_t$aligned_reads <- unlist(lapply(explodeCigarOpLengths(init_t$cigar, ops = c("M", "I")), function(x) sum(x)))
  
  # calc read identity..
  init_t_final <- init_t %>%
    dplyr::mutate(identity = (1 - NM/aligned_reads)*100) %>%
    dplyr::group_by(minion_read_name) %>%
    dplyr::filter(identity == max(identity),
                  aligned_reads == max(aligned_reads)) %>%
    dplyr::distinct(minion_read_name, .keep_all = T) %>%
    dplyr::mutate(gene = str_split_fixed(mapped_gene,"-",2)[,2])
  
  return(init_t_final)
}


group_maps16 <- function(input,bc){
  
  border <- 500
  m1 <-  border - 6 
  m2 <- (1599741)-(hvo_gff_rrna$end_feature[1]-border)+2
  m3 <- (1599741)-(hvo_gff_rrna$end_feature[1]-border)+2+(hvo_gff_rrna$start_feature[1]-1598084)
  c1 <- 20
  
  read_in_bam_large(input) %>%
    ungroup() %>%
    dplyr::filter(identity >= 90, strand == "+") %>%
    dplyr::mutate(group = case_when((start >= m2 & start < m3 & end >= m3 & soft_l >= c1) ~ "(01) 5_extended_16S",
                                    (start >= m2 & start < m3 & end >= m3 & soft_l < c1) ~ "(02) post-16S-bhb",
                                    (start <= m1 & end >= m3) ~ "(03) circ_16S",
                                    (start <= m1 & end >= m2 & end < m3 & soft_r < c1) ~ "(04) open_circ_16S_5",
                                    (end <= m1 & soft_r <= c1 & soft_l < 2000) ~ "(06) mature_16S",
                                    (start >= m3 & soft_l <= c1 & soft_r < 2000) ~ "(06) mature_16S")) %>%
    dplyr::mutate(barcode = bc)
}

group_maps23 <- function(input,bc){
  
  border <- 500
  m1 <- border + 33 
  m2 <- (hvo_gff_rrna$end_feature[2]+72)-(hvo_gff_rrna$end_feature[2]-border)+2
  m3 <- (hvo_gff_rrna$end_feature[2]+72)-(hvo_gff_rrna$end_feature[2]-border)+2+(hvo_gff_rrna$start_feature[2]-11-(hvo_gff_rrna$start_feature[2]-50))
  c1 <- 20
  
  mat_l <- (hvo_gff_rrna$end_feature[2]+33)-(hvo_gff_rrna$start_feature[2]-11)
  read_in_bam_large(input) %>%
    ungroup() %>%
    dplyr::filter(identity >= 90, strand == "+") %>%
    dplyr::mutate(group = case_when((start < m1 & end > m3) ~ "(10) circ_23S",
                                    (start < (m1) & end > m1 & end <= m2 & soft_r <= c1) ~ "(11) open_circ_23S_cutBHB3",
                                    (start < m1  & end <= m1 & soft_r <= c1 & width < mat_l) ~ "(12) mature_23S",
                                    (start >= m3  & end > m3 & soft_l <= c1 & width < mat_l) ~ "(12) mature_23S")) %>%
    dplyr::mutate(barcode = bc)
}

group_maps_all <- function(input,bc){
  
  bhb1 <- 1598084
  bhb2 <- 1599741
  bhb3 <- hvo_gff_rrna$start_feature[2] - 50 - 11
  bhb4 <- hvo_gff_rrna$end_feature[2] + 72
  mat1 <- hvo_gff_rrna$start_feature[1]
  mat2 <- hvo_gff_rrna$end_feature[1]-8
  mat3 <- hvo_gff_rrna$start_feature[2] - 11
  mat4 <- hvo_gff_rrna$end_feature[2] + 32
  mat5 <- hvo_gff_rrna$start_feature[3]
  mat6 <- hvo_gff_rrna$end_feature[3]
  trna_end <- 1599833
  c1 <- 20
  
  read_in_bam_large(input) %>%
    ungroup() %>%
    dplyr::filter(identity >= 90, strand == "+") %>%
    dplyr::mutate(group = case_when((start < bhb1 & start >= (mat1-444) & end > mat1 & soft_l <= c1 & soft_r <= c1) ~ "(01) 5_extended_16S",
                                    (start >= bhb1 & start < mat1 & end > mat1 & end <= bhb2 & soft_l <= c1 & soft_r <= c1) ~ "(02) post-16S-bhb",
                                    (start >= mat1 & end > mat2 & end <= bhb2 & soft_r < c1 & soft_l < c1) ~ "(05) open_circ_16S_cutBHB3",
                                    (start >= mat1 & start < mat2 & end > mat1 & end <= mat2 & soft_l <= c1 & soft_r <= c1) ~ "(06) mature_16S",
                                    (start >= trna_end & start < bhb3 & end > mat3 & soft_l <= c1 & soft_r <= c1) ~ "(07) 5_extended_23S",
                                    (start >= trna_end & end > bhb4 & end <= mat5 & soft_l <= c1 & soft_r <= c1) ~ "(08) 3_extended_23S",
                                    (start >= bhb3 & start < mat3 & end > mat4 & end <= bhb4 & soft_l < c1 & soft_r < c1) ~ "(09) post-23S-bhb",
                                    (start >= mat3 & start < mat4 & end > mat4 & end <= bhb4 & soft_l < c1 & soft_r < c1) ~ "(11) open_circ_23S_cutBHB3",
                                    (start >= mat3 & end <= mat4 & soft_l <= c1 & soft_r <= c1) ~ "(12) mature_23S")) %>%
    dplyr::mutate(barcode = bc)
}

write_fastq_group_lines <- function(set,g){

  index_g <- parse_number(str_split_fixed(g, " ", 2)[,1])

  out_read <- comb_map_summary %>%
    dplyr::filter(group_all %in% g,
                  barcode.x == set |barcode == set | barcode.y == set ) %>%
    dplyr::select(minion_read_name) %>%
    deframe()
  dir.create(paste0(dir, "split_fastq/", set),showWarnings = F)
  write_lines(x = out_read,
              file = paste0(dir, "split_fastq/",set, "/", set, "_",index_g,".lst"))
}
                                    
# data ----
wan_bc <- paste0(rep("barcode0",4),c(1:3,9))

## genome ====
hvo_gff_rrna <- read_in_gff_rrna(here("data/genome/hvo.gff"))

## circ category table ====
dir <- ".path_to_cDNA_analysis"

### 16S ####
list16circ <- list.files(paste0(dir,"mapped_fl_only_circ16S"), recursive = T, pattern = ".sorted.bam$",full.names = T)
list16barc <- str_extract(list16circ, "barcode\\d+")
circ_map16_cat <- pmap_dfr(list(list16circ[c(1:3,9)],(list16barc[c(1:3,9)])),group_maps16)

### 23S ####
list23circ <- list.files(paste0(dir,"mapped_fl_only_circ23S"), recursive = T, pattern = ".sorted.bam$",full.names = T)
list23barc <- str_extract(list23circ, "barcode\\d+")
circ_map23_cat <- pmap_dfr(list(list23circ[1:4],(list23barc[1:4])),group_maps23)

## normal mapping table ====
listallcirc <- list.files(paste0(dir,"mapped_fl_only"), recursive = T, pattern = ".sorted.bam$",full.names = T)
listallbarc <- str_extract(listallcirc, "barcode\\d+")
all_map_cat <- pmap_dfr(list(listallcirc[c(1:3,9)],(listallbarc[c(1:3,9)])),group_maps_all)

## combine tables ====
comb_map <- all_map_cat %>%
  full_join(circ_map16_cat %>%
              dplyr::rename(group_16c = group) %>%
              dplyr::select(barcode,minion_read_name, group_16c),
            by = "minion_read_name") %>%
  full_join(circ_map23_cat %>%
              dplyr::rename(group_23c = group) %>%
              dplyr::select(barcode,minion_read_name, group_23c),
            by = "minion_read_name") %>%
  dplyr::mutate(group_all = ifelse(!is.na(group_16c), group_16c, group_23c),
                group_all = ifelse(is.na(group_all),group, group_all),
                clip_g = case_when(((soft_l > soft_r) & (soft_l > 20)) ~ "left-clipped",
                                   ((soft_l < soft_r) & (soft_l > 20)) ~ "right-clipped",
                                   ((soft_r > soft_l) & (soft_r > 20)) ~ "right-clipped",
                                   ((soft_r < soft_l) & (soft_r > 20)) ~ "left-clipped",
                                   (soft_l <= 20) ~ "none", 
                                   (soft_r <= 20) ~ "none"),
                group_all = ifelse((group == "(05) open_circ_16S_cutBHB3" & group_all == "(06) mature_16S"), "(05) open_circ_16S_cutBHB3",
                                   ifelse((group_all == "(01) 5_extended_16S" & (soft_r > 20 | soft_l > 20)), NA,
                                           ifelse((group_all == "(02) post-16S-bhb" & (soft_r > 20 | soft_l > 20)), NA,
                                                         ifelse((group_all == "(06) mature_16S" & (soft_r > 20 | soft_l > 20)),NA,
                                                                ifelse((group_all == "(11) open_circ_23S_cutBHB3" & (clip_g != "none")),NA,
                                                                       ifelse((group_all == "(12) mature_23S" & (clip_g != "none")),NA,group_all)))))))

## write fastq read groups to output ====
wanted <- c(levels(as.factor(comb_map_summary$group_all))[c(1,2,3,4,5,6)])
for(i in 1:6){write_fastq_group_lines("barcode01", wanted[i])}



# exploratory plots ----
## start-end plot ====
non_circ <-c(levels(as.factor(comb_map$group_all))[c(1,2,5,6,7,8,9,11,12)], NA)


ggplot(data = comb_map %>%
         dplyr::filter(start > hvo_gff_rrna$start_feature[1]-500,
                       end < hvo_gff_rrna$end_feature[2]+500,
                       barcode.x == "barcode01" | barcode.y == "barcode01" | barcode == "barcode01",
                       group_all %in% non_circ) %>%
         dplyr::mutate(clip_g = factor(clip_g, levels = c("none", "left-clipped", "right-clipped"))) %>%
         dplyr::filter(!is.na(clip_g)) %>%
         group_by(group_all) %>%
         arrange(desc(start)) %>%
         dplyr::mutate(rown = 1:n())) +
  facet_grid(rows = vars(group_all), scales = "free_y") +
  geom_segment(aes(y = rown, yend = rown, x = start, xend = end, color = clip_g), size = 1) +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  ylab("") +
  scale_color_manual(values = c("#A2B2C1","#0C365D", "#9B882D")) +
  geom_hline(yintercept = c(0)) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())



## Circ plotting ====
### Circ 16S ####
circ16 <- levels(as.factor(comb_map$group_all))[c(3,4)]

border <- 500
m1 <-  border - 6 
m2 <- (1599741)-(hvo_gff_rrna$end_feature[1]-border)+2
m3 <- (1599741)-(hvo_gff_rrna$end_feature[1]-border)+2+(hvo_gff_rrna$start_feature[1]-1598084)

ggplot(data = circ_map16_cat %>%
         dplyr::filter(barcode == "barcode01",
                       group %in% circ16) %>%
         dplyr::mutate(clip_g = case_when(((soft_l > soft_r) & (soft_l > 20)) ~ "left-clipped",
                                          ((soft_l < soft_r) & (soft_l > 20)) ~ "right-clipped",
                                          ((soft_r > soft_l) & (soft_r > 20)) ~ "right-clipped",
                                          ((soft_r < soft_l) & (soft_r > 20)) ~ "left-clipped",
                                          (soft_l <= 20) ~ "none", 
                                          (soft_r <= 20) ~ "none")) %>%
         dplyr::mutate(clip_g = factor(clip_g, levels = c("none", "left-clipped", "right-clipped"))) %>%
         group_by(group) %>%
         arrange(desc(start)) %>%
         dplyr::mutate(rown = 1:n())) +
  facet_grid(rows = vars(group), scales = "free_y") +
  geom_vline(xintercept = c(m1,m2,m3),
             linetype = "dashed") +
  geom_segment(aes(y = rown, yend = rown, x = start, xend = end, color = clip_g), size = .5) +
  theme_minimal() +
  ylab("") +
  scale_color_manual(values = c("#A2B2C1","#0C365D", "#9B882D")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

### Circ 23S ####
circ23 <- levels(as.factor(comb_map$group_all))[c(10)]
border <- 500
m1_23 <- border + 33 
m2_23 <- (hvo_gff_rrna$end_feature[2]+72)-(hvo_gff_rrna$end_feature[2]-border)+2
m3_23 <- (hvo_gff_rrna$end_feature[2]+72)-(hvo_gff_rrna$end_feature[2]-border)+2+(hvo_gff_rrna$start_feature[2]-11-(hvo_gff_rrna$start_feature[2]-50))

ggplot(data = circ_map23_cat %>%
         dplyr::filter(barcode == "barcode01",
                       group %in% circ23) %>%
         dplyr::mutate(clip_g = case_when(((soft_l > soft_r) & (soft_l > 20)) ~ "left-clipped",
                                          ((soft_l < soft_r) & (soft_l > 20)) ~ "right-clipped",
                                          ((soft_r > soft_l) & (soft_r > 20)) ~ "right-clipped",
                                          ((soft_r < soft_l) & (soft_r > 20)) ~ "left-clipped",
                                          (soft_l <= 20) ~ "none", 
                                          (soft_r <= 20) ~ "none")) %>%
         group_by(group) %>%
         dplyr::mutate(clip_g = factor(clip_g, levels = c("none", "left-clipped", "right-clipped"))) %>%
         arrange(desc(start)) %>%
         dplyr::mutate(rown = 1:n())) +
  facet_grid(rows = vars(group), scales = "free_y") +
  geom_vline(xintercept = c(m1_23,m2_23,m3_23),
             linetype = "dashed") +
  geom_segment(aes(y = rown, yend = rown, x = start, xend = end, color = clip_g), size = .5) +
  theme_minimal() +
  ylab("") +
  scale_color_manual(values = c("#A2B2C1","#0C365D", "#9B882D")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

# Stats ----
## reads mapping to 16S and 23S ====
comb_map %>%
  dplyr::mutate(rRNA_group = ifelse(str_detect(group_all, "23S"), "23S",
                                    ifelse(str_detect(group_all, "16S"), "16S","else"))) %>%
  replace_na(list(rRNA_group = "else")) %>%
  dplyr::mutate(rRNA_group = factor(rRNA_group, levels = rev(c("16S", "23S", "else")))) %>%
  dplyr::mutate(bc = ifelse(!is.na(barcode.x),barcode.x,
                            ifelse(is.na(barcode.x) & !is.na(barcode.y), barcode.y, barcode))) %>%
  dplyr::filter(bc %in% c("barcode01", "barcode02", "barcode09")) %>%
  group_by(bc) %>%
  dplyr::mutate(total_reads = n()) %>%
  group_by(bc, rRNA_group) %>%
  summarise(frac = n()/total_reads*100) %>%
  distinct(frac, bc, rRNA_group) %>%
  ggplot(aes(y = fct_reorder(bc,rev(bc)), x = frac, fill = rRNA_group)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("grey45", "#D5D9E5","#8591B3")) +
  geom_text(aes(label = round(frac,0)),
            position = position_stack(vjust = 0.5), size = 3) +
  scale_x_continuous(expand = c(0,0)) +
  theme_minimal() +
  ylab("") +
  theme(panel.grid.minor = element_blank())
  
## 16S categories ====
range16S <- c((hvo_gff_rrna$start_feature[1]-500):(hvo_gff_rrna$end_feature[1]+200))
map_na_16S <- comb_map %>%
  dplyr::filter(barcode.x %in% paste0(rep("barcode",3), c("01", "02", "09")),
                start %in% range16S, end %in% range16S, is.na(group_all)) %>%
  group_by(barcode.x) %>%
  summarise(counts = n()) %>%
  dplyr::rename(bc = 1) %>%
  dplyr::mutate(group_all = NA)

comb_map %>%
  dplyr::mutate(rRNA_group = ifelse(str_detect(group_all, "23S"), "23S",
                                    ifelse(str_detect(group_all, "16S"), "16S","else"))) %>%
  dplyr::filter(rRNA_group %in% c("16S")) %>%
  dplyr::mutate(bc = ifelse(!is.na(barcode.x),barcode.x,
                            ifelse(is.na(barcode.x) & !is.na(barcode.y), barcode.y, barcode))) %>%
  dplyr::filter(bc %in% c("barcode01", "barcode02", "barcode09")) %>%
  group_by(bc, group_all) %>%
  dplyr::mutate(bc = factor(bc, levels = rev(paste0(rep("barcode",3), c("01", "02", "09"))))) %>%
  summarise(counts = n()) %>%
  distinct(bc, group_all, counts) %>%
  bind_rows(map_na_16S) %>%
  group_by(bc) %>%
  dplyr::mutate(frac = counts/sum(counts)*100) %>%
  ggplot(aes(y = (bc), x = frac, fill = group_all)) +
  geom_bar(stat = "identity", color = "black", size = 0.5) +
  scale_x_continuous(expand = c(0,0)) +
  theme_minimal() +
  ylab("") +
  theme(panel.grid.minor = element_blank())  +
  scale_fill_manual(values = rev(c("#9A9A9A", "#99ffce","#62e3c6","#78a3ec","#7b53df","#5c1d9f")),na.value = "white") 

## 23S categories ====
range23S <- c((hvo_gff_rrna$start_feature[2]-200):(hvo_gff_rrna$end_feature[2]+200))
map_na_23S <- comb_map %>%
  dplyr::filter(barcode.x %in% paste0(rep("barcode",3), c("01", "02")),
                start %in% range23S, end %in% range23S, is.na(group_all)) %>%
  group_by(barcode.x) %>%
  summarise(counts = n()) %>%
  dplyr::rename(bc = 1) %>%
  dplyr::mutate(group_all = NA)

comb_map %>%
  dplyr::mutate(rRNA_group = ifelse(str_detect(group_all, "23S"), "23S",
                                    ifelse(str_detect(group_all, "16S"), "16S","else"))) %>%
  dplyr::filter(rRNA_group %in% c("23S")) %>%
  dplyr::mutate(bc = ifelse(!is.na(barcode.x),barcode.x,
                            ifelse(is.na(barcode.x) & !is.na(barcode.y), barcode.y, barcode))) %>%
  dplyr::filter(bc %in% c("barcode01", "barcode02")) %>%
  group_by(bc, group_all) %>%
  dplyr::mutate(bc = factor(bc, levels = rev(paste0(rep("barcode",3), c("01", "02", "09"))))) %>%
  summarise(counts = n()) %>%
  distinct(bc, group_all, counts) %>%
  bind_rows(map_na_23S) %>%
  group_by(bc) %>%
  dplyr::mutate(frac = counts/sum(counts)*100) %>%
  ggplot(aes(y = (bc), x = frac, fill = group_all)) +
  geom_bar(stat = "identity", color = "black", size = 0.5) +
  scale_x_continuous(expand = c(0,0)) +
  theme_minimal() +
  ylab("") +
  theme(panel.grid.minor = element_blank()) +
  scale_fill_manual(values = rev(c("#CFCFCF", "#916ECA","#D05FAD","#F2698B","#FD9E7E","#EDBA5E")), 
                    na.value = "white") 
```

> Important

Note that for classification of circular precursors (class 3) using
direct RNA reads, the group_maps16 function was slightly modified to
account for lower sequencing quality.
`(start <= m1 & end >= m3) ~ "(03) circ_16S"`, was therefore replaced
with `(start <= m1 & end >= m3+5) ~ "(03) circ_16S"`.

### Circular RNA detection

To investigate circular rRNA reads in more detail, permuted linear
sequences were created. These sequences contained 500 nt upstream of the
annotated rRNA end to the predicted 3´cleavage site of the bhb site and
were joined with the 5´-clevage site of the bhb up to 500 nt downstream
of the annotated rRNA start.

``` r
# Pyrococcus furiosus ====
## pfu genome/fasta information ====
pfu_gff <- read.gff(paste0(dir,"genome/pfu.gff")) %>%
  dplyr::filter(type == "rRNA")
pfu_fasta <- readDNAStringSet(paste0(dir,"genome/pfu.fasta"))
names(pfu_fasta)[1] <- "chr"

## write post 16S bhb ====
cut16_pfu <- paste(pfu_fasta$chr[(pfu_gff$start[1]-500):(120624)],
                   pfu_fasta$chr[122213:(pfu_gff$end[2]+500)],
                          sep = "")

cut16_pfu_all <- paste(pfu_fasta$chr[(pfu_gff$start[1]-500):(120624)],
                       pfu_fasta$chr[122213:122422],
                       pfu_fasta$chr[125518:(125518+300)],
                       sep = "")


### write permuted 16S circ-rRNA sequence ####
pfu_circ16_junc_open <- paste(pfu_fasta$chr[(pfu_gff$end[1]-500):(122213)],
                              pfu_fasta$chr[120624:(pfu_gff$start[1]+500)],
                              sep = "")

### gff type annotation for IGV visualisation ####
gff_circ16S_pfu <- data.table(seqid = "pfu_circ_16S_open",
                              source = "custom",
                              type = "gene",
                              start = c(1,492,522, 573),
                              end = c(491,521,572,1074),
                              score = ".",
                              strand = "+",
                              phase = ".",
                              attributes = c("16S_end", "trailBHB", "leadBHB", "16S_start"))
## 23S ====
# > 120624 5´end BHB
# > 122213 3´end BHB
pfu_gff_rRNA$start_feature[2]-29
pfu_gff_rRNA$end_feature[2]+166

### write permuted 23S circ-rRNA sequence ####
pfu_circ23_junc_open <- paste(pfu_fasta$CP023154[(pfu_gff_rrna$end_feature[2]-500):(125517)],
                               pfu_fasta$CP023154[122422:(pfu_gff_rRNA$start_feature[2]+500)],
                               sep = "")

### gff type annotation for IGV visualisation ####
gff_circ23S_pfu <- data.table(seqid = "pfu_circ_23S_open",
                              source = "custom",
                              type = "gene",
                              start = c(1,499,668, 697),
                              end = c(498,667,696,1198),
                              score = ".",
                              strand = "+",
                              phase = ".",
                              attributes = c("23S_end", "trailBHB", "leadBHB", "23S_start"))

# write to fasta - PFU ----
## circ 16S ====
vroom_write(gff_circ16S_pfu, file = paste0(dir,"genome/circ/circ16S_pfu.gff"), col_names = F, delim = "\t")
write.fasta(sequences = pfu_circ16_junc_open, 
            names = "pfu_circ_16S_open",
            file.out = paste0(dir,"genome/circ/circ16S_pfu.fasta"))

## circ 23S ====
vroom_write(gff_circ23S_pfu, file = paste0(dir,"genome/circ/circ23S_pfu.gff"), col_names = F, delim = "\t")
write.fasta(sequences = pfu_circ23_junc_open, 
            names = "pfu_circ_23S_open",
            file.out = paste0(dir,"genome/circ/circ23S_pfu.fasta"))


# Sulfolobus acidocaldarius ----
## sac genome/fasta information ====
sac_gff <- read.gff(paste0(dir,"genome/sac.gff")) %>%
  dplyr::filter(type == "rRNA")
sac_fasta <- readDNAStringSet(paste0(dir,"genome/sac.fasta"))
names(sac_fasta)[1] <- "chr"

## write post 16S bhb ====
cut16_sac <- paste(reverseComplement(sac_fasta$chr[(1108646):(sac_gff$end[2]+500)]),
                   reverseComplement(sac_fasta$chr[(sac_gff$start[1]-500):(1107096)]),
                          sep = "")

## write post 16S bhb ====
cut16_hvo <- paste(fasta$chr[(gff$start[1]-500):(1598084)],
                   fasta$chr[1599741:(gff$end[2]+500)],
                          sep = "")

### write permuted 16S circ-rRNA sequence ####
sac_circ16_junc_open <- paste(reverseComplement(sac_fasta$chr[(1107096):(sac_gff$start[2]+500)]),reverseComplement(sac_fasta$chr[(sac_gff$end[2]-500):(1108646)]),
                          sep = "")

gff_circ16S_sac <- data.table(seqid = "sac_circ_16S_open",
                              source = "custom",
                              type = "gene",
                              start = c(1,500,543, 559),
                              end = c(499,542,558,1060),
                              score = ".",
                              strand = "+",
                              phase = ".",
                              attributes = c("16S_end", "trailBHB", "leadBHB", "16S_start"))

### write permuted 23S circ-rRNA sequence ####
# > 1106952 5´end BHB
# > 1103876 3´end BHB
mat3 <- sac_gff_rRNA$end_feature[1] - 16
mat4 <- sac_gff_rRNA$start_feature[1] + 8
bhb3 <- 1106952
bhb4 <- 1103876
sac_circ23_junc_open <- paste(reverseComplement(sac_fasta$chr[(bhb4):(sac_gff$start[1]+500)]),
      reverseComplement(sac_fasta$chr[(sac_gff$end[1]-500):(bhb3)]),
                          sep = "")


gff_circ23S_sac <- data.table(seqid = "sac_circ_23S_open",
                              source = "custom",
                              type = "gene",
                              start = c(1,493,523, 556),
                              end = c(492,522,555,1198),
                              score = ".",
                              strand = "+",
                              phase = ".",
                              attributes = c("23S_end", "trailBHB", "leadBHB", "23S_start"))

## post16Sbhb ====
write.fasta(sequences = cut16_sac, 
            names = "cut16_sac",
            file.out = paste0(dir,"genome/cut/cut16S_sac.fasta"))

## circ 16S ====
vroom_write(gff_circ16S_sac, file = paste0(dir,"genome/circ/circ16S_sac.gff"), col_names = F, delim = "\t")
write.fasta(sequences = sac_circ16_junc_open, 
            names = "sac_circ_16S_open",
            file.out = paste0(dir,"genome/circ/circ16S_sac.fasta"))

## circ 23S ====
vroom_write(gff_circ23S_sac, file = paste0(dir,"genome/circ/circ23S_sac.gff"), col_names = F, delim = "\t")
write.fasta(sequences = sac_circ23_junc_open, 
            names = "sac_circ_23S_open",
            file.out = paste0(dir,"genome/circ/circ23S_sac.fasta"))

# Haloferax volcanii -----
## haloferax genome/fasta information ====
gff <- read.gff(paste0(dir,"genome/hvo.gff")) %>%
  dplyr::filter(type == "rRNA")
fasta <- readDNAStringSet(paste0(dir,"genome/hvo.fasta"))
names(fasta)[1] <- "chr"


## write post 16S bhb ====
# > leading16S-23S cut at BHB < #
# > 1598084 5´end BHB
# > 1599741 3´end BHB

## write permuted 16S circ-rRNA sequence ====
circ16_junc_open <- paste(fasta$chr[(gff$end[1]-500):(1599741)],
                          fasta$chr[1598084:(gff$start[1]+500)],
                          sep = "")

gff_circ16S_hvo <- data.table(seqid = "circ_16S_open",
                              source = "custom",
                              type = "gene",
                              start = c(1,502,571, 679),
                              end = c(501,570,678,1179),
                              score = ".",
                              strand = "+",
                              phase = ".",
                              attributes = c("16S_end", "trailBHB", "leadBHB", "16S_start"))

## 23S ==== 
# > 1598084 5´end BHB
# > 1599741 3´end BHB
circ23_junc_open <- paste(fasta$chr[(gff$end[2]-500):(gff$end[2]+72)],
                          fasta$chr[(gff$start[2]-50):(gff$start[2]+500)],
                          sep = "")

gff_circ23S_hvo <- data.table(seqid = "circ_23S_open",
                              source = "custom",
                              type = "gene",
                              start = c(1,534,574,613),
                              end = c(533,573,612,1125),
                              score = ".",
                              strand = "+",
                              phase = ".",
                              attributes = c("23S_end", "trailBHB23", "leadBHB23", "23S_start"))


# write to fasta - HVO ----

## circ 16S ====
vroom_write(gff_circ16S_hvo, file = paste0(dir,"genome/circ/circ16S_hvo.gff"), col_names = F, delim = "\t")
write.fasta(sequences = circ16_junc_open, 
            names = "circ_16S_open",
            file.out = paste0(dir,"genome/circ/circ16S_hvo.fasta"))

## circ 23S ====
vroom_write(gff_circ23S_hvo, file = paste0(dir,"genome/circ/circ23S_hvo.gff"), col_names = F, delim = "\t")
write.fasta(sequences = circ23_junc_open, 
            names = "circ_23S_open",
            file.out = paste0(dir,"genome/circ/circ23S_hvo.fasta"))
```

Nanopore reads were re-mapped to the linear permuted sequences and again
categorised by their 5´ ends and 3´ ends as circular (reads cover
unconnected mature parts of the rRNA and extend over the bhb) or
opened-circular pre-rRNAs (read extends over bhb, 3´break at mature rRNA
start) (*compare detection of rRNA processing*).

### Modified base detection

#### Stage processing

For stage-dependent modified base detection direct RNA reads were sorted
according to the read groups defined using the direct cDNA reads with
two minor modifications. The quality threshold was lowered from 90 % to
80 % and the soft-clipping cutoff was increased from 20 to 30
nucleotides to account for the lower sequencing quality of direct RNA
sequencing using Nanopore technology.

``` r
### R ###

# load libraries ----
library(here)
source(here("Rscripts/load_libraries.R"))

# functions ----
read_in_bam_large <- function(file){
  
  init <- readGAlignments(file, use.names = T, param = ScanBamParam(tag=c("NM"), what="mapq"))
  init_t <- GenomicAlignments::as.data.frame(init) %>%
    dplyr::mutate(minion_read_name = names(init),
                  mapped_gene = seqnames) 
  
  left  <- paste(str_split_fixed(string = init_t$cigar, pattern = "M", n = 2)[,1],"M", sep = "")
  right <- paste(str_split_fixed(string = init_t$cigar, pattern = "M", n = 2)[,2],"1M", sep = "")
  
  #................................calculate cigar tables / SOFT AND HARD CLIPPING!!!
  init_t$soft_l <- as_tibble(cigarOpTable(left))$S
  init_t$hard_l <- as_tibble(cigarOpTable(left))$H
  init_t$soft_r <- as_tibble(cigarOpTable(right))$S
  init_t$hard_r <- as_tibble(cigarOpTable(right))$H
  init_t$matches <- as_tibble(cigarOpTable(init_t$cigar))$M
  init_t$insertions <- as_tibble(cigarOpTable(init_t$cigar))$I
  init_t$deletions <- as_tibble(cigarOpTable(init_t$cigar))$D
  init_t$skip <- as_tibble(cigarOpTable(init_t$cigar))$N
  init_t$cigar_S <- as_tibble(cigarOpTable(init_t$cigar))$S
  
  # calculate number of aligned reads based on CIGAR operations (M,I)
  init_t$aligned_reads <- unlist(lapply(explodeCigarOpLengths(init_t$cigar, ops = c("M", "I")), function(x) sum(x)))
  
  # calc read identity..
  init_t_final <- init_t %>%
    dplyr::mutate(identity = (1 - NM/aligned_reads)*100) %>%
    dplyr::group_by(minion_read_name) %>%
    dplyr::filter(identity == max(identity),
                  aligned_reads == max(aligned_reads)) %>%
    dplyr::distinct(minion_read_name, .keep_all = T) %>%
    dplyr::mutate(gene = str_split_fixed(mapped_gene,"-",2)[,2])
  
  return(init_t_final)
}


group_maps16 <- function(input,bc, myIdentity, myClip ){
  
  border <- 500
  m1 <-  border - 6 
  m2 <- (1599741)-(hvo_gff_rrna$end_feature[1]-border)+2
  m3 <- (1599741)-(hvo_gff_rrna$end_feature[1]-border)+2+(hvo_gff_rrna$start_feature[1]-1598084)
  c1 <- myClip
  
  read_in_bam_large(input) %>%
    ungroup() %>%
    dplyr::filter(identity >= myIdentity, strand == "+") %>%
    dplyr::mutate(group = case_when((start >= m2 & start < m3 & end >= m3 & soft_l >= c1) ~ "(01) 5_extended_16S",
                                    (start >= m2 & start < m3 & end >= m3 & soft_l < c1) ~ "(02) post-16S-bhb",
                                    (start <= m1 & end >= m3+5) ~ "(03) circ_16S",
                                    (start <= m1 & end >= m2 & end < m3 & soft_r < c1) ~ "(04) open_circ_16S_5",
                                    (start <= m1 & end >= m1 & end <= m2 & soft_r < c1) ~ "(05) open_circ_16S_cutBHB3",
                                    (end <= m1 & soft_r <= c1 & soft_l < 2000) ~ "(06) mature_16S",
                                    (start >= m3 & soft_l <= c1 & soft_r < 2000) ~ "(06) mature_16S")) %>%
    dplyr::mutate(barcode = bc)
}

group_maps23 <- function(input,bc, myIdentity, myClip ){
  
  border <- 500
  m1 <- border + 33 
  m2 <- (hvo_gff_rrna$end_feature[2]+72)-(hvo_gff_rrna$end_feature[2]-border)+2
  m3 <- (hvo_gff_rrna$end_feature[2]+72)-(hvo_gff_rrna$end_feature[2]-border)+2+(hvo_gff_rrna$start_feature[2]-11-(hvo_gff_rrna$start_feature[2]-50))
  c1 <- myClip
  mat_l <- (hvo_gff_rrna$end_feature[2]+33)-(hvo_gff_rrna$start_feature[2]-11)
  
  read_in_bam_large(input) %>%
    ungroup() %>%
    dplyr::filter(identity >= myIdentity, strand == "+") %>%
    dplyr::mutate(group = case_when((start < m1 & end > m3) ~ "(10) circ_23S",
                                    (start < (m1) & end > m1 & end <= m2 & soft_r <= c1) ~ "(11) open_circ_23S_cutBHB3",
                                    (start < m1  & end <= m1 & soft_r <= c1 & width < mat_l) ~ "(12) mature_23S",
                                    (start >= m3  & end > m3 & soft_l <= c1 & width < mat_l) ~ "(12) mature_23S")) %>%
    dplyr::mutate(barcode = bc)
}

group_maps_all <- function(input,bc, myIdentity, myClip ){
  
  bhb1 <- 1598084
  bhb2 <- 1599741
  bhb3 <- hvo_gff_rrna$start_feature[2] - 50 - 11
  bhb4 <- hvo_gff_rrna$end_feature[2] + 72
  mat1 <- hvo_gff_rrna$start_feature[1]
  mat2 <- hvo_gff_rrna$end_feature[1]
  mat3 <- hvo_gff_rrna$start_feature[2] - 11
  mat4 <- hvo_gff_rrna$end_feature[2] + 32
  mat5 <- hvo_gff_rrna$start_feature[3]
  mat6 <- hvo_gff_rrna$end_feature[3]
  trna_end <- 1599833
  c1 <- myClip
  
  read_in_bam_large(input) %>%
    ungroup() %>%
    dplyr::filter(identity >= myIdentity, strand == "+") %>%
    dplyr::mutate(group = case_when((start < bhb1 & start >= (mat1-444) & end > mat1 & soft_l <= c1 & soft_r <= c1) ~ "(01) 5_extended_16S",
                                    (start >= bhb1 & start < mat1 & end > mat2 & end <= bhb2 & soft_l <= c1 & soft_r <= c1) ~ "(02) post-16S-bhb",
                                    (start >= mat1 & end > mat2 & end <= bhb2 & soft_r < c1 & soft_l < c1) ~ "(05) open_circ_16S_cutBHB3",
                                    (start >= mat1 & start < mat2 & end > mat1 & end <= mat2 & soft_l <= c1 & soft_r <= c1) ~ "(06) mature_16S",
                                    (start >= trna_end & start < bhb3 & end > mat3 & soft_l <= c1 & soft_r <= c1) ~ "(07) 5_extended_23S",
                                    (start >= trna_end & end > bhb4 & end <= mat5 & soft_l <= c1 & soft_r <= c1) ~ "(08) 3_extended_23S",
                                    (start >= bhb3 & start < mat3 & end > mat4 & end <= bhb4 & soft_l < c1 & soft_r < c1) ~ "(09) post-23S-bhb",
                                    (start >= mat3 & start < mat4 & end > mat4 & end <= bhb4 & soft_l < c1 & soft_r < c1) ~ "(11) open_circ_23S_cutBHB3",
                                    (start >= mat3 & end <= mat4 & soft_l <= c1 & soft_r <= c1) ~ "(12) mature_23S")) %>%
    dplyr::mutate(barcode = bc)
}

write_fastq_group_lines <- function(set,g){
  index_g <- parse_number(str_split_fixed(g, " ", 2)[,1])
  out_read <- comb_mapDRS %>%
    dplyr::filter(group_all %in% g,
                  barcode == set) %>%
    dplyr::select(minion_read_name) %>%
    deframe()
  dir.create(paste0(dir3, "split_fastq/", set),showWarnings = F)
    write_lines(x = out_read,
                file = paste0(dir3, "split_fastq/",set, "/", set, "_",index_g,".lst"))
}

# data ----

## genome data ====
### only rRNAs ####
hvo_gff_rrna   <- read_in_gff_rrna(paste0(dir, "/genome/hvo.gff")) 

### including rRNAs and tRNAs ####
hvo_gff_s <- read.gff(here("data/genome/hvo.gff")) %>%
  dplyr::filter(start > hvo_gff_rrna$start_feature[1]-500,
                end < hvo_gff_rrna$start_feature[3]+200) %>%
  distinct(start, end, .keep_all = T)

## circ category table ====
### 16S ####
dir2 <- "path_to_DRS_analsis"
list16circDRS <- list.files(paste0(dir2,"mapped_fl_only_circ16S"), recursive = T, pattern = ".sorted.bam$",full.names = T)
list16barcDRS <- str_split_fixed(list16circDRS, "/", 8)[,7]
circ_map16_catDRS <- pmap_dfr(list(list16circDRS[c(1:2)],(list16barcDRS[c(1:2)]), 80, 30),group_maps16)

### 23S ####
list23circDRS <- list.files(paste0(dir2,"mapped_fl_only_circ23S"), recursive = T, pattern = ".sorted.bam$",full.names = T)
list23barcDRS <- str_split_fixed(list23circDRS, "/", 8)[,7]
circ_map23_catDRS <- pmap_dfr(list(list23circDRS[c(1:2)],(list23barcDRS[c(1:2)]), 80, 30),group_maps23)

## normal mapping table ====
listallcircDRS <- list.files(paste0(dir2,"mapped_fl_only"), recursive = T, pattern = ".sorted.bam$",full.names = T)
listallbarcDRS <- str_split_fixed(listallcircDRS, "/", 8)[,7]
all_map_catDRS <- pmap_dfr(list(listallcircDRS[c(1:2)],(listallbarcDRS[c(1:2)]), 80, 30),group_maps_all)

## combine tables ====
comb_mapDRS <- all_map_catDRS %>%
  ungroup() %>%
  full_join(circ_map16_catDRS %>%
              ungroup() %>%
              dplyr::rename(group_16c = group) %>%
              dplyr::select(barcode,minion_read_name, group_16c),
            by = "minion_read_name") %>%
  full_join(circ_map23_catDRS %>%
              ungroup() %>%
              dplyr::rename(group_23c = group) %>%
              dplyr::select(barcode,minion_read_name, group_23c),
            by = "minion_read_name") %>%
  dplyr::mutate(group_all = ifelse(!is.na(group_16c), group_16c, group_23c),
                group_all = ifelse(is.na(group_all),group, group_all),
                clip_g = case_when(((soft_l > soft_r) & (soft_l > 30)) ~ "left-clipped",
                                   ((soft_l < soft_r) & (soft_l > 30)) ~ "right-clipped",
                                   ((soft_r > soft_l) & (soft_r > 30)) ~ "right-clipped",
                                   ((soft_r < soft_l) & (soft_r > 30)) ~ "left-clipped",
                                   (soft_l <= 30) ~ "none", 
                                   (soft_r <= 30) ~ "none"),
                group_all = ifelse((group == "(05) open_circ_16S_cutBHB3" & group_all == "(06) mature_16S"), "(05) open_circ_16S_cutBHB3",
                                   ifelse((group_all == "(01) 5_extended_16S" & (soft_r > 30 | soft_l > 30)), NA,
                                          ifelse((group_all == "(02) post-16S-bhb" & (soft_r > 30 | soft_l > 30)), NA,
                                                 ifelse((group_all == "(06) mature_16S" & (soft_r > 30 | soft_l > 30)),NA,
                                                        ifelse((group_all == "(11) open_circ_23S_cutBHB3" & (clip_g != "none")),NA,
                                                               ifelse((group_all == "(12) mature_23S" & (clip_g != "none")),NA,group_all)))))))


## write read groups to id files ====
group_16S <- levels(as.factor(comb_mapDRS$group_all))[grep("16S",levels(as.factor(comb_mapDRS$group_all)))]
for(i in 1:6){write_fastq_group_lines("hvo_notex", group_16S[i])}
for(i in 1:6){write_fastq_group_lines("hvo_notex_dksga", group_16S[i])}


# exploratory ----
## start-end plot ====
levels(as.factor(comb_mapDRS$barcode))

all_list <- levels(as.factor(comb_map$group_all))
S16_list <- all_list[1:6]
S23_list <- c(all_list[7:12])

ggplot(data = comb_mapDRS %>%
         dplyr::filter(start > hvo_gff_rrna$start_feature[1]-500,
                       end < hvo_gff_rrna$end_feature[2]+500,
                       barcode.x == "hvo_notex") %>%
         group_by(group_all) %>%
         arrange(desc(start)) %>%
         dplyr::mutate(rown = 1:n())) +
  geom_vline(xintercept = hvo_gff_s$start) +
  geom_vline(xintercept = hvo_gff_s$end) +
  facet_grid(rows = vars(group_all), scales = "free_y") +
  geom_segment(aes(y = rown, yend = rown, x = start, xend = end, color = clip_g), size = 1) +
  theme_minimal() +
  ylab("") +
  scale_color_manual(values = c("#224156", "#98AAB9", "#F1AC96"))+
  geom_vline(xintercept = c(1598084,1599741), linetype = "dashed") +
  geom_vline(xintercept = c(1600007,1602997), linetype = "dashed") +
  geom_vline(xintercept = c(hvo_gff_rrna$start_feature[1]-444,
                            hvo_gff_rrna$start_feature[1]-287,
                            hvo_gff_rrna$start_feature[1]-185), linetype = "dashed") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) 
  coord_cartesian(xlim =  c(hvo_gff_rrna$end_feature[1]-20,hvo_gff_rrna$end_feature[1]+100))
dev.off()


# New analysis ----

## numbers in group =====
### 16S categories ####
range16S <- c((hvo_gff_rrna$start_feature[1]-500):(hvo_gff_rrna$end_feature[1]+200))
map_na_16S <- comb_mapDRS %>%
  dplyr::filter(barcode.x %in% c("hvo_notex", "hvo_notex_dksga"),
                start %in% range16S, end %in% range16S, is.na(group_all)) %>%
  group_by(barcode.x) %>%
  summarise(counts = n()) %>%
  dplyr::rename(bc = 1) %>%
  dplyr::mutate(group_all = NA)

comb_mapDRS %>%
  dplyr::mutate(rRNA_group = ifelse(str_detect(group_all, "23S"), "23S",
                                    ifelse(str_detect(group_all, "16S"), "16S","else"))) %>%
  dplyr::filter(rRNA_group %in% c("16S")) %>%
  dplyr::mutate(bc = ifelse(!is.na(barcode.x),barcode.x,
                            ifelse(is.na(barcode.x) & !is.na(barcode.y), barcode.y, barcode))) %>%
  dplyr::filter(bc %in% c("hvo_notex", "hvo_notex_dksga")) %>%
  group_by(bc, group_all) %>%
  summarise(counts = n()) %>%
  distinct(bc, group_all, counts) %>%
  bind_rows(map_na_16S) %>%
  group_by(bc) %>%
  dplyr::mutate(frac = counts/sum(counts)*100) %>%
  ggplot(aes(y = (bc), x = frac, fill = group_all)) +
  geom_bar(stat = "identity", color = "black", size = 0.5) +
  scale_x_continuous(expand = c(0,0)) +
  theme_minimal() +
  ylab("") +
  theme(panel.grid.minor = element_blank())  +
  scale_fill_manual(values = rev(c("#9A9A9A", "#99ffce","#62e3c6","#78a3ec","#7b53df","#5c1d9f")),na.value = "white") 

### 23S categories ####
range23S <- c((hvo_gff_rrna$start_feature[2]-200):(hvo_gff_rrna$end_feature[2]+200))
map_na_23S <- comb_mapDRS %>%
  dplyr::filter(barcode.x %in% c("hvo_notex", "hvo_notex_dksga"),
                start %in% range23S, end %in% range23S, is.na(group_all)) %>%
  group_by(barcode.x) %>%
  summarise(counts = n()) %>%
  dplyr::rename(bc = 1) %>%
  dplyr::mutate(group_all = NA)

comb_mapDRS %>%
  dplyr::mutate(rRNA_group = ifelse(str_detect(group_all, "23S"), "23S",
                                    ifelse(str_detect(group_all, "16S"), "16S","else"))) %>%
  dplyr::filter(rRNA_group %in% c("23S")) %>%
  dplyr::mutate(bc = ifelse(!is.na(barcode.x),barcode.x,
                            ifelse(is.na(barcode.x) & !is.na(barcode.y), barcode.y, barcode))) %>%
  dplyr::filter(bc %in% c("hvo_notex", "hvo_notex_dksga")) %>%
  group_by(bc, group_all) %>%
  summarise(counts = n()) %>%
  distinct(bc, group_all, counts) %>%
  bind_rows(map_na_23S) %>%
  group_by(bc) %>%
  dplyr::mutate(frac = counts/sum(counts)*100) %>%
  ggplot(aes(y = (bc), x = frac, fill = group_all)) +
  geom_bar(stat = "identity", color = "black", size = 0.5) +
  scale_x_continuous(expand = c(0,0)) +
  theme_minimal() +
  ylab("") +
  theme(panel.grid.minor = element_blank()) +
  scale_fill_manual(values = rev(c("#CFCFCF", "#916ECA","#D05FAD","#F2698B","#FD9E7E","#EDBA5E")), 
                    na.value = "white") 
```

#### Stage sorting

Reads IDs were extracted for each stage and direct RNA fastq files
sorted using seqtk subseq (v. 1.3-r106).

#### ESB calculation

Sorted Fastq files were subsequently mapped and used for calculating the
frequency of correct, deleted, inserted and wrong nucleotides at each
genomic position using pysamstats (v. v1.1.2).

``` bash
# run pysamstats (polyA-trimmed direct RNA-seq data as input)
for file in ${ext_dir}/analysis/mapped_cut/*/*.sorted.bam
do 
  filename_extended=${file##*/}
  foldername=$(echo $filename_extended | cut -d"_" -f 1)
  ext=$(echo $filename_extended | cut -d"." -f 1)
  echo ${ext} pileup started ...

  out_py=${ext_dir}/analysis/pileups/${ext}
  mkdir -p ${out_py}
  
  pysamstats \
  -D 8000000 \
  -S nofilter \
  --window-size=1 \
  -t variation_strand \
  -f ${ext_dir}/analysis/genome/hvo.fasta \
  $file > ${out_py}/${ext}_pileup.tsv
done
```

> Note that the sum of substitution, deletions and insertion frequencies
> was defined as the Error of Specific Bases (% ESB).

``` r
### R ###
# functions ----
pile_esb_all <- function(input, grouptype){
  pile_base_all(input, grouptype) %>%  
    dplyr::filter(category %in% c("correct", "wrong", "deletion", "insertion")) %>%
    group_by(pos) %>%
    dplyr::mutate(total = sum(n)) %>%
    ungroup() %>%
    dplyr::filter(category != "correct") %>%
    group_by(pos) %>% 
    dplyr::mutate(esb_sum = sum(n)) %>%
    distinct(pos, total, esb_sum) %>%
    dplyr::mutate(frac = esb_sum/total * 100,
                  type = grouptype)
}

# data ----
pileup_files <- list.files(paste0(dir,"pileups/all"), 
                           recursive = T, pattern = "pileup.tsv.gz$", full.names = T)
samples <- data.table(type = 1:4,
                      name = c("wt", "wt", "dksga", "dksga"),
                      rep = c(1,2,1,2))

## calc stats ====
pile_esbs <- pmap_dfr(list(pileup_files_f[c(1:4)],1:4),pile_esb_all) %>%
  left_join(samples) 


### Figure 6 C ####
pil_3c <- pile_esbs %>%
  dplyr::filter(total > 500,
                pos %in% (hvo_gff_rrna$start_feature[1]+13):(hvo_gff_rrna$end_feature[1]-8)) %>%
  group_by(pos, name) %>%
  summarise(mean_frac = mean(frac)) %>%
  pivot_wider(names_from = name, values_from = mean_frac) %>%
  drop_na() %>%
  dplyr::mutate(quot = (wt/dksga)) %>%
  dplyr::filter(!is.na(quot)) %>% 
  ungroup()

# plotting ----
ggplot(data = pil_3c,
       aes(x = pos, y = (quot))) +
  geom_line() +
  geom_area() +
  theme_pubclean() +
  ylab("") +
  geom_vline(xintercept = c(hvo_gff_rrna$start_feature[1]+13,
                            hvo_gff_rrna$end_feature[1]-8),
             linetype = "dashed") +
  geom_hline(yintercept = c(1)) +
  scale_y_continuous(limits = c(0,10))
```

#### Eligos2

Additionally, modified base detection was performed using Eligos2 in
pair_diff_mod (–oddR 0 –esb 0 –pval 1 –adjPval 1) to identify RNA
modifications in the wildtype sample compared to the ∆KsgA sample
(control sample).

First, a custom bed file (16S region) was created in R:

``` r
### R ###
hvo_gff_16S <- ape::read.gff(paste0(dir,"/genome/hvo.gff")) %>%
  dplyr::filter(type %in% "rRNA", str_detect(attributes, "16S")) %>%
  slice(1)

custom_bed <- data.table(chrom  = hvo_gff_16S$seqid,
                         start  = hvo_gff_16S$start,
                         end    = hvo_gff_16S$end,
                         name   = hvo_gff_16S$type,
                         score  = hvo_gff_16S$score,
                         strand = hvo_gff_16S$strand)

vroom_write(x = custom_bed, file = "direct_rna_data/genome/hvo_16S.bed", col_names = F)
```

Next, Eligos2 was performed according to the description at
<https://gitlab.com/piroonj/eligos2>.

``` bash
## Index reference sequence
samtools faidx hvo.fasta

## Run ELIGOS compare between samples when Wild-type (-tbam) and Knock-out (-cbam) for all pairwise comparisons
eligos2 pair_diff_mod \
  -tbam hvo_wt_replicate1.sorted.bam hvo_wt_replicate2.sorted.bam \
  -cbam hvo_dksga_replicate1.sorted.bam hvo_dksga_replicate2.sorted.bam \
  -reg ksga_data/hvo_16S.bed \
  -ref ksga_data/hvo.fasta \
  --oddR 0 \
  --esb 0 \
  --pval 1 \
  --adjPval 1 \
  -t 8 \
  -o ksga_data/results/
```

Plotting was done in R:

``` r
### R ###

# libraries ----
source(here("Rscripts/load_libraries.R"))

# functions ----
read_in_eli <- function(precursor){
  vroom(paste0(ext_dir, precursor, "eligos_output_folder/hvo_notex.sorted_vs_hvo_notex_dksga.sorted_on_hvo_16S_combine.txt")) %>%
    dplyr::mutate(pre = as.factor(precursor))
}

base_annotation <- function(coord_left, coord_right, fasta_file){
  fasta <- readDNAStringSet(fasta_file)
  names(fasta) <- "genome"
  sequence <- as.character(fasta$genome[coord_left:coord_right])
  return(sequence)
}

# data ----
ext_dir <- "/...eligos_output"

## eligos data ====
eligos_f <- list.files(dir,
                      recursive = T, pattern = "_combine.txt", full.names = T)
eligos_samples <- data.table(type = 1:length(eligos_f),
                             precursor = str_remove_all(str_split_fixed(eligos_f, "\\/", 8)[,7], "_new"),
                             rep_wt = str_remove_all(str_split_fixed(str_split_fixed(eligos_f, "\\/", 8)[,8],".sorted",2)[,1], "bc_"),
                             rep_dksga = str_remove_all(str_split_fixed(str_split_fixed(eligos_f, "\\/", 8)[,8],"\\.",4)[,2], "sorted_vs_bc_"))

eligos_frame <- pmap_dfr(list(eligos_f,1:nrow(eligos_samples)), read_in_eli) %>%
  left_join(eligos_samples) 

## genome ====
hvo_gff_rRNA   <- read_in_gff_rrna(paste0(dir, "/genome/hvo.gff")) 
interesting_postions_hvo <- c(1598192+910, 1598192+1352, 1598192+1432, 1598192+1450,1598192+1451, 1598192+1442)

## filter for interesting region ====
eligos_frame_f <- eligos_frame %>%
  dplyr::rename(pos = end_loc) %>% 
  dplyr::select(pos, oddR, adjPval, total_reads, precursor, rep_wt, rep_dksga) 
eligos_frame_f2 <- eligos_frame_f %>%
  dplyr::mutate(oddR = ifelse(oddR == Inf, 0, oddR)) %>%
  group_by(precursor,pos) %>%
  mutate(adjPval_BH = p.adjust(adjPval, method = "BH")) %>%
  summarize(mean_oddR = weighted.mean(oddR,w = total_reads, na.rm = T),
            signif = sum(adjPval_BH < 0.01)/n())

## e.g. H45 ====
set_p <- 20
int_p <- 4
int_seq <- base_annotation(coord_left = (interesting_postions_hvo[int_p]-set_p), 
                           coord_right = (interesting_postions_hvo[int_p]+set_p),
                           fasta_file = here("data/genome/hvo.fasta"))
 
eligos_frame_f3 <-  eligos_frame_f2 %>%
  dplyr::filter(pos > interesting_postions_hvo[int_p]-set_p,
                pos < interesting_postions_hvo[int_p]+set_p) %>%
  group_by(precursor) %>%
  complete(pos = full_seq(pos,1), 
           fill = list(mean_oddR = 0, signif = 0)) %>%
  group_by(pos, precursor) %>%
  dplyr::mutate(lower = 0,
                upper = max(mean_oddR, na.rm = T))

# plotting ----
ggplot(data = eligos_frame_f3,
       aes(x = pos, y = mean_oddR,group = precursor, 
           fill = signif == 1)) +
  facet_grid(rows = vars(precursor)) +
  scale_fill_manual(values = c("#D3D3D3","#B293F4")) +
  geom_ribbon(aes(min = lower, ymax = upper), 
              fill = "grey85") +
  geom_line(size = 1, color = "black") +
  geom_point(shape = 21, color = "black", size = 1.5) +
  coord_cartesian(xlim = c(1599642-13, 1599642+13)) +
  scale_x_continuous(breaks = (1599642-13):(1599642+13),
                     labels = strsplit(as.character(hvo_fasta$chr[(1599642-13):(1599642+13)]),"*")[[1]]) +
  scale_y_continuous(limits = c(0,50)) +
  theme_pubclean() +
  ylab("Odds ratio") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())
dev.off()
```

#### Modified base detection based on signal data

For modified base detection based on signal data, Fast5 files were first
sorted to pre-determined rRNA stages using fast5_subset (ont_fast5_api).
Next, multi-Fast5 files were converted to single read files using
multi_to_single_fast5 (ont_fast5_api). Finally reads were pre-processed,
resquiggled and the raw signals used to extract signal level and dwell
value information using tombo (Version 1.5.1,
<https://nanoporetech.github.io/tombo>). The tombo pipeline was run for
pre-sorted files as well as for unsorted files (Stoiber et al. 2016).
Downstream analysis was performed using custom R scripts (R Foundation
for Statistical Computing. 2018) similar to the ones used for Eligos
analysis.
