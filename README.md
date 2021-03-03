# Hewill Aligner
Hewill Aligner is a highly optimized Burrow-Wheeler Aligner specifically for Illumina HiSeq 2500 Rapid Mode WGS short-read alignment.

## Compiler
- GCC >= 10.2
- Intel Threading Building Blocks (`sudo apt install libtbb-dev`)

## Run
- build executable
```
$ g++-10 main.cpp -o hewill -pthread -ltbb -std=c++20 -O3 biomodern/ssw.cpp -Wno-ignored-attributes
```
- index (only support for uncompressed hs37d5.fa, [download])
```
$ ./hewill index /mnt/fa/hs37d5.fa
```
- align (only support for uncompressed fastq)
- The following command will generate *HG001.1.sam*, *HG001.2.sam*, ... *HG001.X.sam* and *HG001.Y.sam* in */mnt/sam/HG001* folder.
```
// fa_path fq1_path fq2_path sam_prefix sample_name(SM) read_group(RGID) thread_num
$ ./hewill align /mnt/fa/hs37d5.fa /mnt/fq/HG001.1.fq /mnt/fq/HG001.2.fq /mnt/sam/HG001 HG001 1 80
```

## Performance
versus [bwa-mem2]:
<img src="https://raw.githubusercontent.com/hewillk/aligner/master/performance.png" />

## Important
This aligner is highly optimized on the following sequencing characteristic (other datasets are not recommended):
- read length: **~148bp**
- insert size: **~550bp**

[bwa-mem2]: https://github.com/bwa-mem2/bwa-mem2
[download]: https://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
