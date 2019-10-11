longranger align --id=rep2 --fastqs=./ --reference=/public/home/xyhuang/Genome/10Xgenome/dm3/refdata-dm3 --localcores=36

python runMulti.py -i possorted_bam.bam -f summary.csv -s 1 -j /public/home/xyhuang/Tools/littletools/juicer_tools.1.8.9_jcuda.0.8.jar -g ~/Genome/10Xgenome/dm3/dm3.chrom.sizes

perl -lane '@tmp=split(/\t/,$_,3);@tmp2=split(/;/,$tmp[2]);print $_ foreach(@tmp2);' output/out.frag.bed | sort --parallel=4 -k1,1 -k2,2n -k3,3n |bedtools genomecov -bga -i - -g ~/Genome/10Xgenome/dm3/dm3.chrom.sizes > output/out.frag.bedgraph
bedGraphToBigWig output/out.frag.bedgraph ~/Genome/10Xgenome/dm3/dm3.chrom.sizes output/out.frag.bw
computeMatrix scale-regions -S output/out.frag.bw -R ~/Genome/10Xgenome/dm3/dm3.bed -a 2000 -b 2000 -m 2000 --skipZeros --outFileName output/out.frag.profile.gz --numberOfProcessors 4
plotHeatmap -m output/out.frag.profile.gz -out output/out.frag.profile.pdf --colorMap GnBu

# tian：不仅仅考虑了barcode ID，也考虑了测序的机器，*bed.q30.lenBE50.CHR_CLEAN.5e0_3e500.EXT.BCline.GAPEXT3000.FRAGMENT_GEMLIN后续的文件都能反应
# eg：H:fragnumH1|human-10000-100000002-AAACACCAGAAACCTABX1-E00515;chr11      75891303        75891934        2       human-10000-100000002-AAACACCAGAAACCTABX1-E00515-FRAG0000
# eg：H-A:fragnumH2-chrnum1|human-10000-100000107-AAACACCAGCTTACACBX1-E00515;chr2     32916229        32916792        1       human-10000-100000107-AAACACCAGCTTACACBX1-E00515-FRAG0000;chr2  235142786       235145352       13      human-10000-100000107-AAACACCAGCTTACACBX1-E00515-FRAG0001
# 我的只考虑了barcode ID
# 也正因为如此，上传到ncbi的data无法用tian的软件跑，可以分不同的run比对，然后将这个run所有的id加上一个相同的机器ID，比如加 E00515: (分号不能少)

## tian：(只保留了primary hit) 我的保留了secondary hit (STEP0401_F2304.run.p)
# samtools view -h -F 2304 $F|samtools view -hbS - > $file.F2304.bam &

## tian：(没有取唯一比对的结果) 我的增加了提取唯一比对的结果 (AS!=XS) (PX_0403_SH_BCLINE.sh)
## tian：(去掉了部分染色体)  我的全部保留下来了
# cat $F|awk '{if($5>=30) print}' > $file.q30
# cat $file.q30|awk '{if($3-$2>=50)print}' > $file.q30.lenBE50
# cat $file.q30.lenBE50 |grep -v "EBV"|grep -v "_random"|grep -v "chrUn"|grep -v "hs"|grep -v "Het"|grep -v "chrU"|grep -v "chrM"|grep -v "extra" > $file.q30.lenBE50.CHR_CLEAN
# 现在的也增加了black chromsome，可以去掉

## tian：(延伸500bp是 +:5-3 -:3-5) 我的是都是从5-3延伸500bp (PX_0403_EXT_R1_3e500.py)
# word=line.rstrip("\n").split("\t")
# S=word[1]
# E=word[2]
# C=word[0]		
# EXTS=int(S)-int(EXT5)
# EXTE=int(E)+int(EXT3)

## tian：(将跨染色体的结果分成了自连接和染色体间链接) 我的是全部考虑 (STEP040503_CLASS_GEM2EAS.run.p/STEP040504_CLASS_E2AS.run.p/STEP040505_COMBINE_SA.run)

## 最终的日志文件中 Reads_q>=30_&_len>=50bp_&_3'ext500、Reads_q>=30_&_len>=50bp_&_3'ext500_&_MAJOR_CHROM 写反了。

## 