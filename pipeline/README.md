

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
    -   <a
        href="#basecalling-demultiplexing-and-trimming-of-direct-cdna-libraries"
        id="toc-basecalling-demultiplexing-and-trimming-of-direct-cdna-libraries">Basecalling,
        demultiplexing and trimming of direct cDNA libraries</a>
    -   <a href="#demultiplexing-of-basecalled-reads-using-guppy_barcoder"
        id="toc-demultiplexing-of-basecalled-reads-using-guppy_barcoder">Demultiplexing
        of basecalled reads using <code>guppy_barcoder</code></a>
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
folder management is shown in the following:

``` bash
rRNA_maturation/
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
    ├── FASTQ
        ├── normal
        ├── full_length
        ├── cutadapt
        └── cutadapt_SSP
    ├── summary
    ├── barcode
    ├── mapped
        ├── raw
        ├── adapter_trimmed
        └── trimmed
    ├── genome
    ├── salmon
    ├── pychopper
        ├── normal
        └── rescued
    ├── tss
        ├── raw
        └── trimmed
    ├── tts
        ├── raw
        └── trimmed
    ├── bed
    └── coverage_data
        ├── raw
        └── trimmed
```

### Basecalling, demultiplexing and trimming of direct RNA libraries

### Basecalling, demultiplexing and trimming of direct cDNA libraries

After sequencing (and despite live-basecalling) all datasets in the
raw_FAST5 📂 were re-basecalled using `guppy` (v. 6.3.2+bb5453e) in
high-accuracy mode with a q-score cutoff of 7. Next, basecalled files
were demultiplexed using the `guppy_barcoder` command from the `guppy`
suite (available in the [ONT community](https://nanoporetech.com) using
default parameters). Demultiplexed files in FASTQ format were written to
the fastq_pass 📂.

``` bash
# files
input=rRNA_maturation/rebasecallling
output_DRS=microbepore/data/FASTQ/normal/run_id # add run id
output_cDNA=microbepore/data/basecalled/run_id # add run id

# Basecalling of DRS files
guppy_basecaller \
--input_path ${input} \ # input path
--save_path ${output_DRS} \ # output path
-c rna_r9.4.1_70bps_hac.cfg  \ # config file: high accuracy RNA
--calib_detect \ # detect calibration spike-in
--reverse_sequence true \ # reverse since sequenced 3´-->5´
--u_substitution true \ # replace U´s with T´s
--compress_fastq \ # compress output
--fast5_out \ # output FAST5
--recursive \ # look for FAST5 recursively in path
--progress_stats_frequency 60 \ # output progress every minute
--chunks_per_runner 256 \ # options for Mk1C
--gpu_runners_per_device 4 \ # options for Mk1C
--num_callers 1 \ # options for Mk1C
-x auto # options for Mk1C

# Basecalling of cDNA files 
guppy_basecaller \
--input_path ${input} \
--save_path ${output_cDNA} \
-c dna_r9.4.1_450bps_hac.cfg \ # config file: high accuracy cDNA 
--compress_fastq \
--fast5_out \
--recursive \
--progress_stats_frequency 60 \
--chunks_per_runner 256 \
--gpu_runners_per_device 4 \
--num_callers 1 \
-x auto
```

> DRS & cDNA runs require different options.  
> Config file selection based on selected accuracy, flowcell version,
> library preparation kit are listed with
> `guppy_basecaller --print_workflows`

With the selected options `guppy` produces fast5_pass, fast5_fail,
fastq, summary and report files that are written to the FASTQ 📁. FASTQ
are not grouped in pass and fail groups since `--min_qscore` is not
enabled. Multiple FASTQs can be merged using
`cat microbepore/data/basecalled/run_id/*.fastq > microbepore/data/basecalled/run_id/run_id.fastq`.

Sequencing summary files are also written to the FASTQ 📂 and are used
during the quality control of the runs and reads. For better viewing
they can be moved to the summary 📂 using
`mv microbepore/data/FASTQ/run_id/sequencing_summary.txt microbepore/data/summary/run_id.txt`

### Demultiplexing of basecalled reads using `guppy_barcoder`

Next, multiplexed cDNA libraries are demultiplexed in a separate step
using `guppy_barcoder`.

``` bash
# files
input=microbepore/data/basecalled/run_id # add run id
output=microbepore/data/FASTQ/normal/run_id # add run id

# Demultiplexing of (PCR-)cDNA files
guppy_barcoder \
--input_path ${input} \
--save_path ${output} \
--config configuration.cfg \
--barcode_kits SQK-PCB109 \
--progress_stats_frequency 60
```

Multiple FASTQs are written to the FASTQ 📂 and can be merged with
e.g. `cat microbepore/data/FASTQ/run_id/barcode01/*.fastq > microbepore/data/FASTQ/run_id/run_id_barcode01.fastq`.
Barcode summary files are written to the FASTQ 📁 and can be moved to
the barcode 📁.

### Read alignment using [`minimap2`](https://github.com/lh3/minimap2)

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
`gffread microbepore/data/genome/NC_000913.3.gff -g microbepore/data/genome/NC_000913.3.fasta -w microbepore/data/genome/NC_000913.3.transcripts.fasta`.

``` bash
# files
input=microbepore/data/FASTQ/normal # input directory with all merged FASTQ files, 1 for each barcode or single DRS run
fasta=microbepore/data/genome/NC_000913.3.fasta # downloaded from GenBank
transcripts=microbepore/data/genomeNC_000913.3.transcripts.fasta # transcripts file made using gffread

# Mapping & Remapping - loop through all FASTQs
for file in ${input}/*/*.fastq
do
  
  # folder and filenames
  f_ex=${file##*/}
  foldername=$(echo ${f_ex} | cut -d"_" -f 1,2,3) # depending on how you name your files 
  filename=${f_ex%%.*}
  
  # make directories
  mkdir microbepore/data/mapped/raw # direct output to mapped folder for raw reads
  mkdir microbepore/data/mapped/raw/${foldername} # run_id
  output=microbepore/data/mapped/raw/${foldername}/${filename} # run_id/barcode_id
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
input=microbepore/data/FASTQ/normal # input directory with all merged FASTQ files, 1 for each barcode or single DRS run

# perform pychopper for all cDNA and (PCR)-cDNA files
for file in ${input}/*/*.fastq
do 

  # folder and filenames
  f_ex=${file##*/}
  foldername=$(echo $f_ex | cut -d"_" -f 1,2,3)
  filename=${f_ex%%.*}
  
  # make directories
  mkdir microbepore/data/pychopper/normal
  mkdir microbepore/data/pychopper/normal/${foldername}
  output=microbepore/data/pychopper/normal/${foldername}/${filename}
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
input=microbepore/data/pychopper/normal # input directory with all merged FASTQ files, 1 for each barcode or single DRS run

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
  output=microbepore/data/pychopper/rescued/${foldername}/${filename}
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
input=microbepore/data/pychopper/

# merge all full-length and rescued reads as full-length
for file in ${input}/normal/*/*/*full_length_output.fastq # both normal and rescued folders
do 
  filename_extended=${file##*/}
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=$(echo $filename_extended | cut -d"_" -f 1,2,3,4,5)

  keyword=$(echo $foldername | cut -d"_" -f 2) # get libary kit ID
  
  mkdir microbepore/data/FASTQ/full_length
  mkdir microbepore/data/FASTQ/full_length/${foldername}
  output=microbepore/data/FASTQ/full_length/${foldername}/${filename}
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
moved to the microbepore/data/FASTQ/full_length folder and adding
\*\_full_length_all\* to the filename.

#### Remove polyA-tails using [`cutadapt`](https://cutadapt.readthedocs.io/en/stable/)

To evaluate the influence of different trimming approaches on the
accuracy of transcript boundary analysis, we applied additional 5´ and
3´ trimming steps using `cutadapt` v3.2.  
To this end, polyA sequences were removed from the 3´ends:

``` bash
# files
input=microbepore/data/FASTQ/full_length # input directory with all merged FASTQ files, 1 for each barcode or single DRS run

for file in ${input}/*/*/*_full_length_all.fastq
do 

  # folder and filenames
  filename_extended=${file##*/}
  keyword=$(echo $filename_extended | cut -d"." -f 2)
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}
  
  mkdir microbepore/data/FASTQ/cutadapt
  mkdir microbepore/data/FASTQ/cutadapt/${foldername}
  output=microbepore/data/FASTQ/cutadapt/${foldername}/${filename}
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
input=microbepore/data/FASTQ/cutadapt

# >  SSP adapter
for file in ${input}/*/*/*cutadapt.fastq
do 
  filename_extended=${file##*/}
  keyword=$(echo $filename_extended | cut -d"." -f 2)
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}
  
  mkdir microbepore/data/FASTQ/cutadapt_SSP
  mkdir microbepore/data/FASTQ/cutadapt_SSP/${foldername}
  output=microbepore/data/FASTQ/cutadapt_SSP/${foldername}/${filename}
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
input=microbepore/data/FASTQ/cutadapt_SSP
fasta=microbepore/data/genome/NC_000913.3.fasta # downloaded from GenBank

# map (pychopper) > polyA_trimmed > SSP trimmed fastqs
for file in ${input}/*/*/*fastq
do 
  filename_extended=${file##*/}
  foldername=$(echo ${filename_extended} | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}

  mkdir microbepore/data/mapped/adapter_trimmed
  mkdir microbepore/data/mapped/adapter_trimmed/${foldername}
  output=microbepore/data/mapped/adapter_trimmed/${foldername}/${filename}
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
input=microbepore/data/mapped/adapter_trimmed
fasta=microbepore/data/genome/NC_000913.3.fasta # downloaded from GenBank
transcripts=microbepore/data/genomeNC_000913.3.transcripts.fasta # transcripts file made using gffread

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
    
    mkdir microbepore/data/mapped/trimmed
    mkdir microbepore/data/mapped/trimmed/${foldername}
    output=microbepore/data/mapped/trimmed/${foldername}/${filename}
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
input=microbepore/data/mapped

# perform tss detection for pychopper auto > cutadapt_polyA > SSP-cutadapt > clipped  or for raw mapped reads
for file in ${input}/trimmed/*/*/*clipped.sorted.bam # ||  for file in ${input}/raw/*/*/*.sorted.bam
do 
  # file and folder names
  filename_extended=${file##*/}
  keyword=$(echo $filename_extended | cut -d"." -f 2)
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}
  
  # make directories
  mkdir microbepore/data/tss/trimmed
  mkdir microbepore/data/tss/trimmed/${foldername}
  output=microbepore/data/tss/trimmed/${foldername}/${filename}
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
input=microbepore/data/mapped

# perform tts detection for pychopper auto > cutadapt_polyA > SSP-cutadapt > clipped  or for raw mapped reads
for file in ${input}/trimmed/*/*/*clipped.sorted.bam # ||  for file in ${input}/raw/*/*/*.sorted.bam
do 
  filename_extended=${file##*/}
  keyword=$(echo $filename_extended | cut -d"." -f 2)
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}
    
  echo ${filename}

  mkdir microbepore/data/tts/trimmed
  mkdir microbepore/data/tts/trimmed
  mkdir microbepore/data/tts/trimmed/${foldername}
  output=microbepore/data/tts/trimmed/${foldername}/${filename}
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
input=microbepore/data/mapped

# calculate coverage over transcripts with TSS and TTS | for pychopper auto > cutadapt > clipped or RAW 
for file in ${input}/trimmed/*/*/*clipped.sorted.bam # ||  for file in ${input}/raw/*/*/*.sorted.bam
do 
  filename_extended=${file##*/}
  keyword=$(echo $filename_extended | cut -d"." -f 2)
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}

  # mk dirs
  mkdir microbepore/data/coverage/trimmed
  mkdir microbepore/data/coverage/trimmed/${foldername}
  output=microbepore/data/coverage/trimmed/${foldername}/${filename}
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