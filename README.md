multi-chia-tools
========
a python program to handles the multi-chia/chia-drop data

REQUIREMENTS
============
python3<br/>
    - pysam<br/>
    - pybedtools<br/>
R<br/>

USAGE of multi-chia-tools
============
    $ python runMulti.py -i possorted_bam.bam -f summary.csv -j juicer_tools.1.8.9_jcuda.0.8.jar -g Genome/dm3/dm3.chrom.sizes

Release Notes
============

Version 1.0.1
--------------------------
1. add convert to .hic format
2. add drop chromosome in black list
3. add just use fq1 data to generate result
4. add just use intra data to generate interaction

Version 1.0.0
--------------------------
1. Implemented core algorithms of process the chiadrop/multi-chia data
