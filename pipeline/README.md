

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
    -   <a
        href="#detection-of-rrna-processing-sites-and-classification-of-rrna-intermediates"
        id="toc-detection-of-rrna-processing-sites-and-classification-of-rrna-intermediates">Detection
        of rRNA processing sites and classification of rRNA intermediates</a>
    -   <a href="#circular-rna-detection"
        id="toc-circular-rna-detection">Circular RNA detection</a>
    -   <a href="#modified-base-detection"
        id="toc-modified-base-detection">Modified base detection</a>
        -   <a href="#identification-of-full-length-reads-using-pychopper"
            id="toc-identification-of-full-length-reads-using-pychopper">Identification
            of full-length reads using <span><code>pychopper</code></span></a>
        -   <a href="#remove-polya-tails-using-cutadapt"
            id="toc-remove-polya-tails-using-cutadapt">Remove polyA-tails using
            <span><code>cutadapt</code></span></a>
        -   <a href="#remove-remaining-ssp-adapter-using-cutadapt"
            id="toc-remove-remaining-ssp-adapter-using-cutadapt">Remove remaining
            SSP adapter using <span><code>cutadapt</code></span></a>
        -   <a href="#mapping-of-trimmed-reads-removing-clips-using-samclip"
            id="toc-mapping-of-trimmed-reads-removing-clips-using-samclip">Mapping
            of trimmed reads, removing clips using <code>samclip</code></a>
    -   <a href="#detection-of-transcript-boundaries"
        id="toc-detection-of-transcript-boundaries">Detection of transcript
        boundaries</a>
        -   <a href="#5end-detection" id="toc-5end-detection">5´end detection</a>
        -   <a href="#3end-detection" id="toc-3end-detection">3´end detection</a>
    -   <a href="#gene-body-coverage-analysis"
        id="toc-gene-body-coverage-analysis">Gene body coverage analysis</a>

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
    ├── pychopper_edlib
    ├── pychopper_edlib_rescue
    
    
    ├── mapped
        ├── raw
        ├── adapter_trimmed
        └── trimmed
    ├── genome
    ├── bed
    └── coverage_data
        ├── raw
        └── trimmed
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
rebasecallling 📂. Multiple FASTQs can be merged using
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

Multiple FASTQs are written to the demultiplexing 📁 and can be merged
with
e.g. `cat rRNA_maturation/demultiplexing/pass/barcode01/*.fastq > rRNA_maturation/analysis/fastq_pass/barcode01.fastq`.
Barcode summary files are written to the rRNA_maturation/demultiplexing
📁 and can be moved to the rRNA_maturation/analysis/summary 📁.

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

Files were mapped to the reference genome from *Escherichia coli* K-12
MG1655 ([GenBank](https://www.ncbi.nlm.nih.gov/nuccore/545778205):
U00096.3) using `minimap2` (Release 2.18-r1015).  
Genome FASTA and GFF3 files have been downloaded from
[GenBank](https://www.ncbi.nlm.nih.gov/nuccore/545778205). Output
alignments in the SAM format were generated with `-ax splice -k14` for
**Nanopore cDNA-seq** and `-ax splice, -uf, -k14` for **DRS** with i)
`-p 0.99`, to return primary and secondary mappings and ii) with `--MD`,
to include the MD tag for calculating mapping identities. Alignment
files were further converted to BAM files, sorted and indexed using
\[`SAMtools`(<https://github.com/samtools/>).  
To analyse single reads in more detail with respect to the RNA type
(mRNA, rRNA, other ncRNA, unspecified) they map to, BAM files were first
converted back to FASTQ using
[`bedtools`](https://bedtools.readthedocs.io/en/latest/) v2.29.2. Next
FASTQ files were remapped to a transcriptome file using `minimap2` with
the previously mentioned parameters to assign single read names with
feature IDs. The transcript file was made using
[`gffread`](https://github.com/gpertea/gffread) with
`gffread rRNA_maturation/data/genome/NC_000913.3.gff -g rRNA_maturation/data/genome/NC_000913.3.fasta -w rRNA_maturation/data/genome/NC_000913.3.transcripts.fasta`.

``` bash
# files
input=rRNA_maturation/data/FASTQ/normal # input directory with all merged FASTQ files, 1 for each barcode or single DRS run
fasta=rRNA_maturation/data/genome/NC_000913.3.fasta # downloaded from GenBank
transcripts=rRNA_maturation/data/genomeNC_000913.3.transcripts.fasta # transcripts file made using gffread

# Mapping & Remapping - loop through all FASTQs
for file in ${input}/*/*.fastq
do
  
  # folder and filenames
  f_ex=${file##*/}
  foldername=$(echo ${f_ex} | cut -d"_" -f 1,2,3) # depending on how you name your files 
  filename=${f_ex%%.*}
  
  # make directories
  mkdir rRNA_maturation/data/mapped/raw # direct output to mapped folder for raw reads
  mkdir rRNA_maturation/data/mapped/raw/${foldername} # run_id
  output=rRNA_maturation/data/mapped/raw/${foldername}/${filename} # run_id/barcode_id
  mkdir ${output}

  if [[ $filename =~ "RNA" ]]; 
  then
  # align using minimap2
  minimap2 -ax splice -p 0.99 -uf -k14 --MD -t 8 ${fasta} ${file} > ${output}/${filename}.sam # DRS
  else
    minimap2 -ax splice -p 0.99 -k14 --MD -t 8 ${fasta} ${file} > ${output}/${filename}.sam # (PCR-)cDNA
  fi
 
  # convert to sorted.bam file
  samtools view -bS ${output}/${filename}.sam -o ${output}/${filename}.bam
  samtools sort ${output}/${filename}.bam -o ${output}/${filename}.sorted.bam
  samtools index ${output}/${filename}.sorted.bam
  
  # bam to fastq for remapping of mapped reads
  bedtools bamtofastq -i ${output}/${filename}.sorted.bam -fq ${output}/${filename}.remapped.fastq
  
  # map again
  if [[ $filename =~ "RNA" ]]; 
  then
  minimap2 -ax splice -p 0.99 -uf -k14 --MD -t 8 ${transcripts} ${output}/${filename}.remapped.fastq > ${output}/${filename}.remapped.sam
  else
    minimap2 -ax splice -p 0.99 -k14 --MD -t 8 ${transcripts} ${output}/${filename}.remapped.fastq > ${output}/${filename}.remapped.sam
 fi
 
  # convert to sorted.bam file
  samtools view -bS ${output}/${filename}.remapped.sam -o ${output}/${filename}.remapped.bam
  samtools sort ${output}/${filename}.remapped.bam -o ${output}/${filename}.remapped.sorted.bam
  samtools index ${output}/${filename}.remapped.sorted.bam
done
```

### Detection of rRNA processing sites and classification of rRNA intermediates

### Circular RNA detection

### Modified base detection

#### Identification of full-length reads using [`pychopper`](https://github.com/nanoporetech/pychopper)

Full-length cDNA reads containing SSP and VNP primers in the correct
orientation were identified using `pychopper` (v.2.5.0) with standard
parameters using the default pHMM backend and autotuned cutoff
parameters estimated from subsampled data. Save output in pychopper 📁.

``` bash
# files
input=rRNA_maturation/data/FASTQ/normal # input directory with all merged FASTQ files, 1 for each barcode or single DRS run

# perform pychopper for all cDNA and (PCR)-cDNA files
for file in ${input}/*/*.fastq
do 

  # folder and filenames
  f_ex=${file##*/}
  foldername=$(echo $f_ex | cut -d"_" -f 1,2,3)
  filename=${f_ex%%.*}
  
  # make directories
  mkdir rRNA_maturation/data/pychopper/normal
  mkdir rRNA_maturation/data/pychopper/normal/${foldername}
  output=rRNA_maturation/data/pychopper/normal/${foldername}/${filename}
  mkdir ${output}

  # perform pychopper using precomputed q
  cdna_classifier.py \
  -r ${output}/${filename}_report.pdf \
  -t 8 \
  -u ${output}/${filename}_unclassified.fastq \
  -w ${output}/${filename}_rescued.fastq \
  -S ${output}/${filename}_stats.txt \
  $file \
  ${output}/${filename}_full_length_output.fastq
done
```

After a first round, a second round of `pychopper` was applied to the
unclassified direct cDNA reads with DCS-specific read rescue enabled.

``` bash
# files
input=rRNA_maturation/data/pychopper/normal # input directory with all merged FASTQ files, 1 for each barcode or single DRS run

# perform pychopper using the -x rescue option for DCS files
for file in ${input}/*unclassified.fastq # only use unclassified reads from first round as input
do 

  # folder and filenames
  filename_extended=${file##*/}
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}
  
  # make directories
  mkdir ${dir}/data/pychopper/rescued
  mkdir ${dir}/data/pychopper/rescued/${foldername}
  output=rRNA_maturation/data/pychopper/rescued/${foldername}/${filename}
  mkdir ${output}
  
  # perfrom pychopper using -X option for native cDNA datasets
  cdna_classifier.py \
  -r ${output}/${filename}_report.pdf \
  -t 8 \
  -x rescue \
  -u ${output}/${filename}_unclassified.fastq \
  -w ${output}/${filename}_rescued.fastq \
  -S ${output}/${filename}_stats.txt \
  $file \
  ${output}/${filename}_full_length_output.fastq
done
```

Reads from rescued and normal folders were merged and used for
subsequent steps.

``` bash
# files
input=rRNA_maturation/data/pychopper/

# merge all full-length and rescued reads as full-length
for file in ${input}/normal/*/*/*full_length_output.fastq # both normal and rescued folders
do 
  filename_extended=${file##*/}
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=$(echo $filename_extended | cut -d"_" -f 1,2,3,4,5)

  keyword=$(echo $foldername | cut -d"_" -f 2) # get libary kit ID
  
  mkdir rRNA_maturation/data/FASTQ/full_length
  mkdir rRNA_maturation/data/FASTQ/full_length/${foldername}
  output=rRNA_maturation/data/FASTQ/full_length/${foldername}/${filename}
  mkdir ${output}
  
  if [[ $keyword =~ "PCB109" ]]; then
    cat $file ${input}/normal/${foldername}/${filename}/${filename}_rescued.fastq > ${output}/${filename}_full_length_all.fastq
  elif [[ $keyword =~ "DCS109" ]]; then
    cat $file ${input}/normal/${foldername}/${filename}/${filename}_rescued.fastq
    ${input}/rescued/${foldername}/${filename}_unclassified/${filename}_unclassified_full_length_output.fastq
    ${input}/rescued/${foldername}/${filename}_unclassified/${filename}_unclassified_rescued.fastq > ${output}/${filename}_full_length_all.fastq
  fi
done
```

For easier handling in the subsequent steps, DRS FASTQ files are also
moved to the rRNA_maturation/data/FASTQ/full_length folder and adding
\*\_full_length_all\* to the filename.

#### Remove polyA-tails using [`cutadapt`](https://cutadapt.readthedocs.io/en/stable/)

To evaluate the influence of different trimming approaches on the
accuracy of transcript boundary analysis, we applied additional 5´ and
3´ trimming steps using `cutadapt` v3.2.  
To this end, polyA sequences were removed from the 3´ends:

``` bash
# files
input=rRNA_maturation/data/FASTQ/full_length # input directory with all merged FASTQ files, 1 for each barcode or single DRS run

for file in ${input}/*/*/*_full_length_all.fastq
do 

  # folder and filenames
  filename_extended=${file##*/}
  keyword=$(echo $filename_extended | cut -d"." -f 2)
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}
  
  mkdir rRNA_maturation/data/FASTQ/cutadapt
  mkdir rRNA_maturation/data/FASTQ/cutadapt/${foldername}
  output=rRNA_maturation/data/FASTQ/cutadapt/${foldername}/${filename}
  mkdir ${output}
  
  # cutadapt
  cutadapt \
    -a "A{10}" \ # trim polyAs longer than 10 bases from the 3´end
    -e 1 \ # allowed error rate
    -j 0 \ # auto-detect cores
    -o ${output}/${filename}.cutadapt.fastq \
    ${file}
done
```

#### Remove remaining SSP adapter using [`cutadapt`](https://cutadapt.readthedocs.io/en/stable/)

Remove remaining SSP sequences from the 5´ends of the cDNA reads using:

``` bash
input=rRNA_maturation/data/FASTQ/cutadapt

# >  SSP adapter
for file in ${input}/*/*/*cutadapt.fastq
do 
  filename_extended=${file##*/}
  keyword=$(echo $filename_extended | cut -d"." -f 2)
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}
  
  mkdir rRNA_maturation/data/FASTQ/cutadapt_SSP
  mkdir rRNA_maturation/data/FASTQ/cutadapt_SSP/${foldername}
  output=rRNA_maturation/data/FASTQ/cutadapt_SSP/${foldername}/${filename}
  mkdir ${output}

  cutadapt \
    -g "TTTCTGTTGGTGCTGATATTGCTGGG" \
    -e 1 \
    -j 0 \
    -o ${output}/${filename}.cutadapt_SSP.fastq \
    ${file}
done
```

#### Mapping of trimmed reads, removing clips using `samclip`

Finally, trimmed reads were mapped using `minimap2` as described before.
Reads with more than 10 clipped bases on either side were removed from
the alignments using [`samclip`](https://github.com/tseemann/samclip)
(v.0.4.0).

1.  Step: Align

``` bash
input=rRNA_maturation/data/FASTQ/cutadapt_SSP
fasta=rRNA_maturation/data/genome/NC_000913.3.fasta # downloaded from GenBank

# map (pychopper) > polyA_trimmed > SSP trimmed fastqs
for file in ${input}/*/*/*fastq
do 
  filename_extended=${file##*/}
  foldername=$(echo ${filename_extended} | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}

  mkdir rRNA_maturation/data/mapped/adapter_trimmed
  mkdir rRNA_maturation/data/mapped/adapter_trimmed/${foldername}
  output=rRNA_maturation/data/mapped/adapter_trimmed/${foldername}/${filename}
  mkdir ${output}

  ## align using minimap2
  if [[ $filename =~ "RNA" ]]; 
  then
  # align using minimap2
  minimap2 -ax splice -p 0.99 -uf -k14 --MD -t 8 ${fasta} ${file} > ${output}/${filename}.sam
  else
    minimap2 -ax splice -p 0.99 -k14 --MD -t 8 ${fasta} ${file} > ${output}/${filename}.sam
  fi
done
```

2.  Step: Remove clipping \> 10 bases

``` bash
input=rRNA_maturation/data/mapped/adapter_trimmed
fasta=rRNA_maturation/data/genome/NC_000913.3.fasta # downloaded from GenBank
transcripts=rRNA_maturation/data/genomeNC_000913.3.transcripts.fasta # transcripts file made using gffread

# remove reads with more than 10 bases that are clipped on either side. 
for file in ${input}/*/*/*.sam
do 
  filename_extended=${file##*/}
  keyword=$(echo $filename_extended | cut -d"." -f 2)
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}
  
  if [[ $keyword =~ "sam" ]]; then
    echo ${foldername}
    echo ${filename}
    echo ${keyword}
    
    mkdir rRNA_maturation/data/mapped/trimmed
    mkdir rRNA_maturation/data/mapped/trimmed/${foldername}
    output=rRNA_maturation/data/mapped/trimmed/${foldername}/${filename}
    mkdir ${output}
  
    # remove mapped reads with a Maximum clip length to allow (10, 5 is default)
    samclip --max 10 --ref ${fasta} < ${file} > ${output}/${filename}.clipped.sam
    
    # convert to sorted.bam file
    samtools flagstat ${output}/${filename}.clipped.sam > ${output}/${filename}.clipped.stats.txt
    samtools view -bS ${output}/${filename}.clipped.sam -o ${output}/${filename}.clipped.bam
    samtools sort ${output}/${filename}.clipped.bam -o ${output}/${filename}.clipped.sorted.bam
    samtools index ${output}/${filename}.clipped.sorted.bam
    
    ## remap fastq converted reads
  bedtools bamtofastq -i ${output}/${filename}.clipped.sorted.bam -fq ${output}/${filename}.remapped.fastq
  
  ## map again
  if [[ $filename =~ "RNA" ]]; 
  then
  # align using minimap2
  minimap2 -ax splice -p 0.99 -uf -k14 --MD -t 8 ${transcripts} ${file} > ${output}/${filename}.remapped.sam
  else
    minimap2 -ax splice -p 0.99 -k14 --MD -t 8 ${transcripts} ${file} > ${output}/${filename}.remapped.sam
  fi
  
  # convert to sorted.bam file
  samtools view -bS ${output}/${filename}.remapped.sam -o ${output}/${filename}.remapped.bam
  samtools sort ${output}/${filename}.remapped.bam -o ${output}/${filename}.remapped.sorted.bam
  samtools index ${output}/${filename}.remapped.sorted.bam
  fi
done
```

### Detection of transcript boundaries

The determination of enriched 5´and 3´ends was carried out in the same
way, but independently of each other, and is briefly explained in the
following: First, strand-specific read ends in bedgraph format were
created from BAM files using
[`bedtools genomecov`](https://bedtools.readthedocs.io/en/latest/) (-5
or -3 option, -bga). Next, the previously published
[`Termseq_peaks`](https://pypi.org/project/termseq-peaks/) script was
used to call peaks for each sample individually without including
replicates (<https://github.com/NICHD-BSPC/termseq-peaks>). This script
is based on `scipy.signal.find_peaks`, which is running in the
background of `Termseq_peaks` with lenient parameters
(prominence=(None,None), width=(1,None), rel_height=0.75). However, we
deliberately used `Termseq_peaks` since its ability to include
replicates by applying an Irreproducible Discovery Rate method which can
be applied to future studies. For end detection, only the leniently
called peaks in the narrowPeak file were used after adding the number of
counts for each position using `bedtools intersect`.

#### 5´end detection

5´end peak calling was performed in the following way:

``` bash
input=rRNA_maturation/data/mapped

# perform tss detection for pychopper auto > cutadapt_polyA > SSP-cutadapt > clipped  or for raw mapped reads
for file in ${input}/trimmed/*/*/*clipped.sorted.bam # ||  for file in ${input}/raw/*/*/*.sorted.bam
do 
  # file and folder names
  filename_extended=${file##*/}
  keyword=$(echo $filename_extended | cut -d"." -f 2)
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}
  
  # make directories
  mkdir rRNA_maturation/data/tss/trimmed
  mkdir rRNA_maturation/data/tss/trimmed/${foldername}
  output=rRNA_maturation/data/tss/trimmed/${foldername}/${filename}
  mkdir ${output}

  # step 1: calculate 5´positions for plus and minus strand
  bedtools genomecov \
    -ibam ${file} \
    -bga \
    -5 \
    -strand + > ${output}/${filename}.plus.bedgraph
  
  bedtools genomecov \
    -ibam ${file} \
    -bga \
    -5 \
    -strand - > ${output}/${filename}.minus.bedgraph
    
  # step 2: termseq peaks
  termseq_peaks ${output}/${filename}.plus.bedgraph ${output}/${filename}.plus.bedgraph --peaks ${output}/${filename}.plus.peaks --strand +
  termseq_peaks ${output}/${filename}.minus.bedgraph ${output}/${filename}.minus.bedgraph --peaks ${output}/${filename}.minus.peaks --strand -
    
  # step 3: add coverage information
  bedtools intersect \
    -wao \
    -a ${output}/${filename}.plus.peaks.oracle.narrowPeak \
    -b ${output}/${filename}.plus.bedgraph \
    > ${output}/${filename}.plus.peaks.oracle.narrowPeak.counts
    
  bedtools intersect \
    -wao \
    -a ${output}/${filename}.minus.peaks.oracle.narrowPeak \
    -b ${output}/${filename}.minus.bedgraph \
    > ${output}/${filename}.minus.peaks.oracle.narrowPeak.counts

done
```

#### 3´end detection

3´end peak calling was performed in the following way:

``` bash
input=rRNA_maturation/data/mapped

# perform tts detection for pychopper auto > cutadapt_polyA > SSP-cutadapt > clipped  or for raw mapped reads
for file in ${input}/trimmed/*/*/*clipped.sorted.bam # ||  for file in ${input}/raw/*/*/*.sorted.bam
do 
  filename_extended=${file##*/}
  keyword=$(echo $filename_extended | cut -d"." -f 2)
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}
    
  echo ${filename}

  mkdir rRNA_maturation/data/tts/trimmed
  mkdir rRNA_maturation/data/tts/trimmed
  mkdir rRNA_maturation/data/tts/trimmed/${foldername}
  output=rRNA_maturation/data/tts/trimmed/${foldername}/${filename}
  mkdir ${output}

  # step 1: calculate 3´positions for plus and minus strand
  bedtools genomecov \
    -ibam ${file} \
    -bga \
    -3 \
    -strand + > ${output}/${filename}.plus.bedgraph
  
  bedtools genomecov \
    -ibam ${file} \
    -bga \
    -3 \
    -strand - > ${output}/${filename}.minus.bedgraph
    
  # step 2: termseq peaks
  termseq_peaks ${output}/${filename}.plus.bedgraph ${output}/${filename}.plus.bedgraph --peaks ${output}/${filename}.plus.peaks --strand +
    termseq_peaks ${output}/${filename}.minus.bedgraph ${output}/${filename}.minus.bedgraph --peaks ${output}/${filename}.minus.peaks --strand -
    
  # step 3: add coverage information
  bedtools intersect \
    -wao \
    -a ${output}/${filename}.plus.peaks.oracle.narrowPeak \
    -b ${output}/${filename}.plus.bedgraph \
    > ${output}/${filename}.plus.peaks.oracle.narrowPeak.counts
    
  bedtools intersect \
    -wao \
    -a ${output}/${filename}.minus.peaks.oracle.narrowPeak \
    -b ${output}/${filename}.minus.bedgraph \
    > ${output}/${filename}.minus.peaks.oracle.narrowPeak.counts

done
```

### Gene body coverage analysis

To assess the impact of trimmings on gene body coverage, a coverage
meta-analysis was performed. First, a transcript file was created for
all genes with an ONT-annotated primary 5´ and 3´ end (see previous
section). Based on this, strand-specific coverage files were created
from the BAM files and coverage analysis performed using a custom R
script.

``` bash
input=rRNA_maturation/data/mapped

# calculate coverage over transcripts with TSS and TTS | for pychopper auto > cutadapt > clipped or RAW 
for file in ${input}/trimmed/*/*/*clipped.sorted.bam # ||  for file in ${input}/raw/*/*/*.sorted.bam
do 
  filename_extended=${file##*/}
  keyword=$(echo $filename_extended | cut -d"." -f 2)
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}

  # mk dirs
  mkdir rRNA_maturation/data/coverage/trimmed
  mkdir rRNA_maturation/data/coverage/trimmed/${foldername}
  output=rRNA_maturation/data/coverage/trimmed/${foldername}/${filename}
  mkdir ${output}

  # calc coverage
  samtools view -F 16 -o temp.sorted.bam ${file} 
  bedtools coverage \
  -d \
  -a ${dir}/data/bed/transcripts.plus.bedgraph \ # bed file of genes with annotated 5´and 3´end
  -b temp.sorted.bam \
  > ${output}/${filename}.plus.coverage
  
  samtools view -f 16 -o temp.sorted.bam ${file} 
  bedtools coverage \
  -d \
  -a ${dir}/data/bed/transcripts.minus.bedgraph \ # bed file of genes with annotated 5´and 3´end
  -b temp.sorted.bam \
  > ${output}/${filename}.minus.coverage
done
```
