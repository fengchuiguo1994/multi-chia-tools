multi-chia-tools
========
a python program to handles the multi-chia/chia-drop data

REQUIREMENTS
============
python3<br/>
    - pysam<br/>
    - pybedtools<br/>
    - scipy<br/>
    - pandas<br/>
    - numpy<br/>
    - rpy2<br/>
R<br/>

USAGE of multi-chia-tools
============
    $ python runMulti.py -i possorted_bam.bam -f summary.csv -j juicer_tools.1.8.9_jcuda.0.8.jar -g Genome/dm3/dm3.chrom.sizes

the result
--------------------------
- $prefix_log.txt           it contaion all static infomation
- $prefix.bed               bed6 format,the 4th column is barcode-reads ID,every reads has extends 500bp in 5'end(from 3' to 5'),chromosome    start   end barcode-reads ID    .   strand
- $prefix.cv.bed            $prefix.bed file change format,Convenient for downstream analysis.barcode-chromosome  start   end read ID
- $prefix.BC.txt            barcode    reads_count
- $prefix.merge.bed         merge the reads to fragment,barcode-chromosome start   end merge_reads_count   merge_reads
- $prefix.all.frag.bed      barcode fragment_count  "chr start end;...;..."
- $prefix.samechrom.frag.bed the file is same as $prefix.all.frag.bed,but it just contain intra-chromosome,the inter-chromosome is split to intra and self.
- $prefix.IGV.bed           the file can import in IGV and show the barcode fragment,but the result looks bad,so I do a new one by R.the file can be transform to bedgraph/bigwig file,and compare to chipseq data.
- $prefix.distance.txt      two fragment distance in a same barcode
- $prefix.distance.pdf      the distance distribution,up is complex==N,down is complex>=N
- $prefix.total.*.bedpe     不同方法计算得到的包含了染色体间交互的4DN格式的结果
- $prefix.intra.*.bedpe     不同方法计算得到的不包含了染色体间交互的4DN格式的结果
- $prefix.total.*.bedpe.hic 4DN格式转成.hic格式
- $prefix.intra.*.bedpe.hic 4DN格式转成.hic格式
- $prefix.chiapet.cluster   类似于chiapet的cluster文件
- $prefix.chiapet.cluster.withpvalue  计算了cluster的显著度
- $prefix.chiapet.cluster.withpvalue.FDR.txt   过滤筛选高可信的cluster

the distance.r result<br/>
![](https://github.com/fengchuiguo1994/multi-chia-tools/blob/2.0.0/img/distance.png)<br/>
the plot_multi_chia.py result<br/>
![](https://github.com/fengchuiguo1994/multi-chia-tools/blob/2.0.0/img/plot_multi_chia.png)<br/>

Release Notes
============

Version 2.0.0
--------------------------
1. add convert to chiapet 
2. add hypergeom to detect Significant cluster(loop)

Version 1.0.1
--------------------------
1. add convert to .hic format
2. add drop chromosome in black list
3. add just use fq1 data to generate result
4. add just use intra data to generate interaction

Version 1.0.0
--------------------------
1. Implemented core algorithms of process the chiadrop/multi-chia data
