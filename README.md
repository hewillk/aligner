# aligner
## Run
```
$ g++-10 main.cpp -o hewill -pthread -ltbb -std=c++20 -O3 biomodern/ssw.cpp -Wno-ignored-attributes
$ ./hewill index /mnt/fa/hs37d5.fa
// fa fq1 fq2 sam_prefix sample_name(SM) read_group(RGID) threads
$ ./hewill align /mnt/fa/hs37d5.fa /mnt/fq/HG001.1.fq /mnt/fq/HG001.2.fq /mnt/sams/HG001 HG001 1 80
```