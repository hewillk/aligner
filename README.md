# Hewill Aligner
Hewill Aligner is a highly optimized Burrow-Wheeler Aligner (implement using modern C++20) specifically for paired-end short read alignment with read length ~150 bp.
**The elapsed time of mapping sequencing data with 50× coverage (2×180GB) is less than 2 hours under 80 cores and ~18 GB memory usage and has comparable performance to [bwa-mem2].**

## Compiler
- GCC >= 10.2
- Intel Threading Building Blocks (`sudo apt install libtbb-dev`)

## Run
### build executable
```
$ g++-10 main.cpp -o hewill -pthread -ltbb -std=c++20 -O3 biomodern/ssw.cpp -Wno-ignored-attributes
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
// fa_path fq1_path fq2_path sam_prefix sample_name(SM) read_group_id(RGID) thread_num
$ ./hewill align /mnt/fa/hs37d5.fa /mnt/fq/HG001.1.fq /mnt/fq/HG001.2.fq /mnt/sam/HG001 HG001 1 80
```

## Performance
PrecisionFDA [Truth Challenge] benchmark versus [bwa-mem2]:
<img src="https://raw.githubusercontent.com/hewillk/aligner/master/performance.png" />

## Important
This aligner is highly optimized on the following sequencing characteristic (other datasets are not recommended):
- read length: **~150bp**

[bwa-mem2]: https://github.com/bwa-mem2/bwa-mem2
[Truth Challenge]: https://precision.fda.gov/challenges/truth
[download]: https://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
