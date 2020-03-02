
# Retrieving Projects {#sra}

Many of your projects will involve working with publically available data from published research. The data itself is typicaly uploaded to a data repository like NCBI or EMBL/EBI. The files are deposited as Sequence Read Archive (SRA) files and are often searchable as bioprojects. Accessing these files can be accomplished using Bioconductor in R.


```r
library(stringr)
library(SRAdb)
library(Biostrings)
```

First, download the database file to your working director via


```r
getSRAdbFile(destdir=getwd(),destfile='SRAmetadb.sqlite.gz',method)
```

Then, we'll make a connection with the database:




```r
sqlfile <- '~/SRAmetadb.sqlite'
sra_con <- dbConnect(SQLite(),sqlfile)
```

Let's assume we want to download the data from the Gevers IBD study that was deposited with the following accession: PRJNA237362. We can see this project on NCBI here: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA237362 . If we click the link next to 'SRA Experiments' and then click one of the sample links on the subsequent page, we'll end up on a page that looks like this: https://www.ncbi.nlm.nih.gov/sra/SRX1418176[accn] . Here' we can access the SRA project code: SRP040765. This is what we need.

Before we continue, let me go over the SRA file types:

* SRA - Accession information that contains the 5 files below
* SRP - Project information and metadata
* SRS - Sample metadata
* SRX - Experiment metadata including library, platform selection, and processing parametes involved in a particular sequencing experiment
* SRR - Sequencing run information
* SRX - Sequence analysis BAM file information

Now that we have the SRP, let's acquire the files we need, specifically the SRA files:


```r
rs <- listSRAfile(c('SRP040765'), sra_con, fileType = 'sra')
str(rs)
```

```
## 'data.frame':	3408 obs. of  5 variables:
##  $ study     : chr  "SRP040765" "SRP040765" "SRP040765" "SRP040765" ...
##  $ sample    : chr  "SRS587325" "SRS587957" "SRS587956" "SRS695027" ...
##  $ experiment: chr  "SRX691644" "SRX692390" "SRX509287" "SRX693364" ...
##  $ run       : chr  "SRR1564534" "SRR1565230" "SRR1215330" "SRR1566410" ...
##  $ ftp       : chr  "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR156/SRR1564534/SRR1564534.sra" "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR156/SRR1565230/SRR1565230.sra" "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR121/SRR1215330/SRR1215330.sra" "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR156/SRR1566410/SRR1566410.sra" ...
```

**rs** is a dataframe containing the SRP, SRS, SRX, and SRR IDs for a given sequencing run, as well as the ftp link to the actual .sra file containing the sequencing information. These are the links we can use to download the entire set of data we need to perform some analysis.

Now, we could export these links and then just iterate through, maybe using bash with wget, to download all of these files. Alternatively, we can do the following:

We'll get a run ID for a run we'ld like to download:


```r
run <- rs$run[1]
run
```

```
## [1] "SRR1564534"
```


If we want the specific SRR (run) information, we do:


```r
run_info <- getSRA(search_terms='SRP040765', out_types=c('run'),sra_con)
str(run_info)
```

```
## 'data.frame':	3408 obs. of  11 variables:
##  $ run_alias      : chr  "A21YN121012.1.Illumina_P7-Fexexoja.screened.bam" "A21YN121012.1.Illumina_P7-Wodejora.screened.bam" "A2WP7130403.1.Illumina_P7-Xakokoxe.screened.bam" "A1UEN121219.1.Illumina_P7-Birarane.screened.bam" ...
##  $ run            : chr  "SRR1564534" "SRR1565230" "SRR1215330" "SRR1566410" ...
##  $ run_date       : chr  "2012-10-12" "2012-10-12" "2013-04-03" "2012-12-19" ...
##  $ updated_date   : chr  "2014-09-04" "2014-09-05" "2019-12-11" "2014-09-06" ...
##  $ spots          : num  120 449 2102 9532 676152 ...
##  $ bases          : num  4.20e+04 1.57e+05 7.36e+05 3.34e+06 1.37e+08 ...
##  $ run_center     : chr  "BI" "BI" "BI" "BI" ...
##  $ experiment_name: logi  NA NA NA NA NA NA ...
##  $ run_url_link   : logi  NA NA NA NA NA NA ...
##  $ run_entrez_link: logi  NA NA NA NA NA NA ...
##  $ run_attribute  : chr  "analysis_type: AssemblyWithoutReference || flowcell_barcode: A21YN || gssr_id: 247611.0 || instrument_name: SL-"| __truncated__ "analysis_type: AssemblyWithoutReference || flowcell_barcode: A21YN || gssr_id: 247629.0 || instrument_name: SL-"| __truncated__ "analysis_type: AssemblyWithoutReference || flowcell_barcode: A2WP7 || gssr_id: 247671.0 || instrument_name: SL-"| __truncated__ "analysis_type: AssemblyWithoutReference || data_type: 16S || flowcell_barcode: A1UEN || gssr_id: 276389.0 || in"| __truncated__ ...
```

and for the SRS (sample) information:


```r
sample_info <- getSRA(search_terms='SRP040765', out_types=c('sample'),sra_con)
str(sample_info)
```

```
## 'data.frame':	1572 obs. of  10 variables:
##  $ sample_alias      : chr  "SKBTI-0107" "SKBTI-0128" "SKBTI-0172" "SKBTI-0594" ...
##  $ sample            : chr  "SRS587325" "SRS587957" "SRS587956" "SRS695027" ...
##  $ taxon_id          : int  408170 408170 408170 408170 408170 408170 408170 408170 408170 408170 ...
##  $ common_name       : chr  NA NA NA NA ...
##  $ anonymized_name   : chr  NA NA NA NA ...
##  $ individual_name   : chr  NA NA NA NA ...
##  $ description       : logi  NA NA NA NA NA NA ...
##  $ sample_url_link   : logi  NA NA NA NA NA NA ...
##  $ sample_entrez_link: logi  NA NA NA NA NA NA ...
##  $ sample_attribute  : chr  "strain: SKBTI-0107 || collection date: missing || geographic location (country and/or sea, region): USA || spec"| __truncated__ "strain: SKBTI-0128 || collection date: missing || geographic location (country and/or sea, region): USA || spec"| __truncated__ "strain: SKBTI-0172 || collection date: missing || geographic location (country and/or sea, region): USA || spec"| __truncated__ "strain: SKBTI-0594 || collection date: missing || geographic location (country and/or sea, region): USA || spec"| __truncated__ ...
```

and SRX (experiment) information:


```r
experiment_info <- getSRA(search_terms='SRP040765', out_types=c('experiment'),sra_con)
str(experiment_info)
```

```
## 'data.frame':	2708 obs. of  27 variables:
##  $ experiment_alias             : chr  "2949006.WR32770.Solexa-122962.A21YN121012.P" "2949006.WR32770.Solexa-122980.A21YN121012.P" "2949006.WR32770.Solexa-123022.A2WP7130403.P" "2949006.WR33991.Solexa-133729.A1UEN121219.P" ...
##  $ experiment                   : chr  "SRX691644" "SRX692390" "SRX509287" "SRX693364" ...
##  $ experiment_title             : chr  "Illumina amplicon sequencing of metagenomic paired-end library 'Solexa-122962' containing sample 'SKBTI-0107'" "Illumina amplicon sequencing of metagenomic paired-end library 'Solexa-122980' containing sample 'SKBTI-0128'" "Illumina amplicon sequencing of metagenomic paired-end library 'Solexa-123022' containing sample 'SKBTI-0172'" "Illumina amplicon sequencing of metagenomic paired-end library 'Solexa-133729' containing sample 'SKBTI-0594'" ...
##  $ study_name                   : logi  NA NA NA NA NA NA ...
##  $ sample_name                  : logi  NA NA NA NA NA NA ...
##  $ design_description           : chr  "Illumina sequencing of human gut metagenome via polymerase chain reaction" "Illumina sequencing of human gut metagenome via polymerase chain reaction" "Illumina sequencing of human gut metagenome via polymerase chain reaction" "Illumina sequencing of human gut metagenome via polymerase chain reaction" ...
##  $ library_name                 : chr  "Solexa-122962" "Solexa-122980" "Solexa-123022" "Solexa-133729" ...
##  $ library_strategy             : chr  "AMPLICON" "AMPLICON" "AMPLICON" "AMPLICON" ...
##  $ library_source               : chr  "METAGENOMIC" "METAGENOMIC" "METAGENOMIC" "METAGENOMIC" ...
##  $ library_selection            : chr  "PCR" "PCR" "PCR" "PCR" ...
##  $ library_layout               : chr  "PAIRED - NOMINAL_SDEV: 0.0E0; NOMINAL_LENGTH: 390; " "PAIRED - NOMINAL_SDEV: 0.0E0; NOMINAL_LENGTH: 390; " "PAIRED - NOMINAL_SDEV: 0.0E0; NOMINAL_LENGTH: 393; " "PAIRED - NOMINAL_SDEV: 0.0E0; NOMINAL_LENGTH: 382; " ...
##  $ library_construction_protocol: logi  NA NA NA NA NA NA ...
##  $ adapter_spec                 : logi  NA NA NA NA NA NA ...
##  $ read_spec                    : chr  "READ_INDEX: 0; READ_LABEL: forward; READ_CLASS: Application Read; READ_TYPE: Forward; BASE_COORD: 1 || READ_IND"| __truncated__ "READ_INDEX: 0; READ_LABEL: forward; READ_CLASS: Application Read; READ_TYPE: Forward; BASE_COORD: 1 || READ_IND"| __truncated__ "READ_INDEX: 0; READ_LABEL: forward; READ_CLASS: Application Read; READ_TYPE: Forward; BASE_COORD: 1 || READ_IND"| __truncated__ "READ_INDEX: 0; READ_LABEL: forward; READ_CLASS: Application Read; READ_TYPE: Forward; BASE_COORD: 1 || READ_IND"| __truncated__ ...
##  $ platform                     : chr  "ILLUMINA" "ILLUMINA" "ILLUMINA" "ILLUMINA" ...
##  $ instrument_model             : chr  "Illumina MiSeq" "Illumina MiSeq" "Illumina MiSeq" "Illumina MiSeq" ...
##  $ instrument_name              : logi  NA NA NA NA NA NA ...
##  $ platform_parameters          : chr  "INSTRUMENT_MODEL: Illumina MiSeq" "INSTRUMENT_MODEL: Illumina MiSeq" "INSTRUMENT_MODEL: Illumina MiSeq" "INSTRUMENT_MODEL: Illumina MiSeq" ...
##  $ sequence_space               : logi  NA NA NA NA NA NA ...
##  $ base_caller                  : logi  NA NA NA NA NA NA ...
##  $ quality_scorer               : logi  NA NA NA NA NA NA ...
##  $ number_of_levels             : logi  NA NA NA NA NA NA ...
##  $ multiplier                   : logi  NA NA NA NA NA NA ...
##  $ qtype                        : logi  NA NA NA NA NA NA ...
##  $ experiment_url_link          : logi  NA NA NA NA NA NA ...
##  $ experiment_entrez_link       : logi  NA NA NA NA NA NA ...
##  $ experiment_attribute         : chr  "analysis_type: AssemblyWithoutReference || gssr_id: 247611.0 || library_type: 16S || lsid: broadinstitute.org:b"| __truncated__ "analysis_type: AssemblyWithoutReference || gssr_id: 247629.0 || library_type: 16S || lsid: broadinstitute.org:b"| __truncated__ "analysis_type: AssemblyWithoutReference || gssr_id: 247671.0 || library_type: 16S || lsid: broadinstitute.org:b"| __truncated__ "analysis_type: AssemblyWithoutReference || data_type: 16S || gssr_id: 276389.0 || library_type: 16S || lsid: br"| __truncated__ ...
```

Using these commands, you should be able to download the .sra files you need along with all corresponding metadata to do analysis. Still, you might be wondering how you get the .fasta files from the .sra file. Well, the easiest way is to use something called the sra toolkit, which can be found here: https://www.ncbi.nlm.nih.gov/books/NBK158900/ . 

So let's say we aimed to extract the sequences from the following sra files: SRR1635768 and SRR1566401. First, we'd download the files:


```r
sra_dir <- tempdir()

sra_fns <- c("SRR1634425","SRR1634428")
for (sra in sra_fns) getSRAfile(sra, sra_con, fileType = 'sra',destDir=sra_dir)
```

```
## Files are saved to: 
## '/tmp/RtmpB5N0L1'
## 
## Files are saved to: 
## '/tmp/RtmpB5N0L1'
```

Then, we'd use the sra toolkit to extract the sequences. Assuming you have downloaded and installed it, we can do the following


```r
sra_output <- tempdir()

sra_files <- list.files(sra_dir,full.names=TRUE,pattern='\\.sra$')

for (i in seq_along(sra_files)) system2('fastq-dump',args=c(sra_files[i],
                                                              '-O', sra_output,
                                                              '--gzip',
                                                              '--clip',
                                                              '--skip-technical',
                                                              '--dumpbase'))
```

And now we can check:


```r
fqs <- list.files(sra_output,full.names=TRUE)

FASTQ <- readDNAStringSet(fqs,format='fastq')
FASTQ
```

Note the header names; they're simply the SRR ids. You'd have to the sra metadata (see above) to match these to specific samples. For example,


```r
sample_info <- getSRA(search_terms='SRR1635768', out_types=c('sample'),sra_con)
str(sample_info)
```

```
## 'data.frame':	1 obs. of  10 variables:
##  $ sample_alias      : chr "SKBTI-0325"
##  $ sample            : chr "SRS734393"
##  $ taxon_id          : int 408170
##  $ common_name       : logi NA
##  $ anonymized_name   : logi NA
##  $ individual_name   : logi NA
##  $ description       : logi NA
##  $ sample_url_link   : logi NA
##  $ sample_entrez_link: logi NA
##  $ sample_attribute  : chr "strain: SKBTI-0325 || collection date: missing || geographic location (country and/or sea, region): USA || spec"| __truncated__
```

### Fastq Dump for Paired End Reads

It should be noted that simply using the fastq-dump command as above works only if your data *does not consist of paired end reads*. If you happen to have paired end reads, then the following arguments must be added to the fastq-dump call:


```r
for (i in seq_along(sra_files)) system2('fastq-dump',args=c(sra_files[i],
                                                              '-O', sra_output,
                                                              '--gzip',
                                                              '--clip',
                                                              '--skip-technical',
                                                              '--dumpbase',
                                                              '--split-files',
                                                              '--readids'))
```

For more information, see https://edwards.sdsu.edu/research/fastq-dump/
