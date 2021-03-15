# Hewill Aligner
Hewill Aligner is a highly optimized Burrow-Wheeler Aligner (implement using modern C++20) specifically for paired-end short read alignment with read length ~148 bp.
**The elapsed time of mapping sequencing data with 50× coverage (2×180GB) is less than 2 hours under 80 cores and ~18 GB memory usage and has comparable performance to [bwa-mem2].**

## Compiler
- GCC >= 10.2
- Intel Threading Building Blocks (`sudo apt install libtbb-dev`)

## Run
### build executable
```
$ g++-10 main.cpp -o hewill -pthread -ltbb -std=c++20 -O3 biomodern/ssw.cpp -Wno-ignored-attributes
```
### Command
```
Command:
         ./hewill index <hs37d5.fa>
         ./hewill align <hs37d5.fa> <in1.fq> <in2.fq> <sam_prefix> <sample_name> <read_group_id> [insert_mean] [insert_var] [thread_num]
Required Arguments:
       <hs37d5.fa>          Reference sequence hs37d5.fa file path
       <in1.fq>             Read1 FASTQ path
       <in2.fq>             Read2 FASTQ path
       <sam_prefix>         SAM file output prefix
       <sample_name>        SAM file SM tag value
       <read_group_id>      SAM file RGID tag value
       
Optional Arguments:
       [insert_mean]        Mean insert size of the paired-end data, default value is 550
       [insert_var]         Variance insert size of the paired-end data, default value is 150
       [thread_num]         Number of threads, default value is return value of std::thread::hardware_concurrency()
```
### index
- Only support for uncompressed hs37d5.fa. ([download])
- **The suffix array sorting time is less than 3 minutes under 80 cores and ~20 GB memory usage.**
```
$ ./hewill index /mnt/fa/hs37d5.fa
```
### align
- Only support for uncompressed fastq.
- The following command will generate *HG001.1.sam*, *HG001.2.sam*, ... *HG001.X.sam* and *HG001.Y.sam* in */mnt/sam/HG001* folder.
```cpp
// <hs37d5.fa> <in1.fq> <in2.fq> <sam_prefix> <sample_name> <read_group_id> [insert_mean] [insert_var] [thread_num]
$ ./hewill align /mnt/fa/hs37d5.fa /mnt/fq/HG001.1.fq /mnt/fq/HG001.2.fq /mnt/sam/HG001 HG001 1 550 150 80
```

## Performance
PrecisionFDA [Truth Challenge] benchmark versus [bwa-mem2]:
<img src="https://raw.githubusercontent.com/hewillk/aligner/master/performance.png" />

## Important
This aligner is highly optimized on the following sequencing characteristic (other datasets are not recommended):
- read length: **148~150bp**

[bwa-mem2]: https://github.com/bwa-mem2/bwa-mem2
[Truth Challenge]: https://precision.fda.gov/challenges/truth
[download]: https://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
