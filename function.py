import pysam
import time
import re
import os
import subprocess
import random
from operator import itemgetter
from bx.intervals.intersection import Intersecter, Interval

'''
function : give runMulti the deal function
create : huang xingyu
version : 1.0.1
begin : 2019-07-11
end : 
'''

def readfile(fin):
    if fin.endswith(".bam"):
        samfile=pysam.AlignmentFile(fin,'rb')
    else:
        samfile=pysam.AlignmentFile(fin)
    return samfile

def bam2bed(samfile,output,prefix,mapq,readlen,minlength,black):
    if not os.path.exists(output):
        os.makedirs(output)
    outbedfile = os.path.join(output,prefix+".bed")
    outbedfile2 = os.path.join(output,prefix+".cv.bed")
    bed = open(outbedfile,'w')
    bed2 = open(outbedfile2,'w')
    
    typeBC = {}
    R1total = R2total = 0
    R1barcode = R2barcode = 0
    R1countmap = R2countmap = 0
    R1uniqmap = R2uniqmap = 0
    R1Q30map = R2Q30map = 0
    R1lenmap = R2lenmap = 0
    R1black = R2black = 0
    total_reads = 0
    nn = 0
    be = time.time()
    
    for line in samfile:
        nn += 1
        if nn%1000000==0:
            print("deal with {0} M reads,use time {1:.5f} s".format(int(nn/1000000),time.time()-be))
            be = time.time()

        if not line.is_secondary:
            total_reads += 1
            # continue
        if line.is_read1:
            R1total += 1
        else:
            R2total += 1
        if line.has_tag('BX'):
            barcode = re.search("(\w){16}",line.get_tag('BX'))[0]
            if barcode not in typeBC:
                typeBC[barcode] = 0
            typeBC[barcode] += 1
            if line.is_read1:
                R1barcode += 1
            else:
                R2barcode += 1
            if not line.is_unmapped:
                if line.is_read1:
                    R1countmap += 1
                else:
                    R2countmap += 1
                if line.get_tag('AS') != line.get_tag('XS'):
                    if line.is_read1:
                        R1uniqmap += 1
                    else:
                        R2uniqmap += 1
                    if line.mapping_quality >= mapq:
                        if line.is_read1:
                            R1Q30map += 1
                        else:
                            R2Q30map += 1
                        if line.reference_length >= readlen:
                            if line.is_read1:
                                R1lenmap += 1
                            else:
                                R2lenmap += 1
                            if line.reference_name not in black:
                                if line.is_read1:
                                    R1black += 1
                                else:
                                    R2black += 1
                                if line.is_reverse:
                                    if line.reference_start <= minlength:
                                        bed.write(line.reference_name+"\t0\t"+str(line.reference_end)+"\t"+barcode+"-"+line.query_name+"\t.\t-\n")
                                        bed2.write(barcode+"-"+line.reference_name+"\t0\t"+str(line.reference_end)+"\t"+line.query_name+"\n")
                                    else:
                                        bed.write(line.reference_name+"\t"+str(line.reference_start-minlength)+"\t"+str(line.reference_end)+"\t"+barcode+"-"+line.query_name+"\t.\t-\n")
                                        bed2.write(barcode+"-"+line.reference_name+"\t"+str(line.reference_start-minlength)+"\t"+str(line.reference_end)+"\t"+line.query_name+"\n")
                                else:
                                    bed.write(line.reference_name+"\t"+str(line.reference_start)+"\t"+str(line.reference_end+minlength)+"\t"+barcode+"-"+line.query_name+"\t.\t+\n")
                                    bed2.write(barcode+"-"+line.reference_name+"\t"+str(line.reference_start)+"\t"+str(line.reference_end+minlength)+"\t"+line.query_name+"\n")
    bed.close()
    bed2.close()
    ## readsID总数      reads条数(包括secondary)    barcode字典     R1total     R2total     R1有barcode数目     R2有barcode数目     R1有 map    R2有 map    R1有 mapuniq    R2有 mapuniq    R1有 map Q30    R2有 map Q30    R1有 map Q30 len    R2有 map Q30 len    R1有 map Q30 len chrom  R2有 map Q30 len chrom
    ## total_reads      nn                          typeBC          R1total    R2total   R1barcode           R2barcode        R1countmap  R2countmap    R1uniqmap        R2uniqmap    R1Q30map         R2Q30map        R1lenmap            R2lenmap             R1black                 R2black
    return(total_reads,nn,typeBC,R1total,R2total,R1barcode,R2barcode,R1countmap,R2countmap,R1uniqmap,R2uniqmap,R1Q30map,R2Q30map,R1lenmap,R2lenmap,R1black,R2black)

def tobed12(mylist,barcode,out):
    if len(mylist) == 1:
        out.write("{0}\t{1}\t{2}\t{3}\t0\t+\t{4}\t{5}\t0\t1\t{6}\t0,\n".format(mylist[0][0],mylist[0][1],mylist[0][2],barcode,mylist[0][1],mylist[0][2],(int(mylist[0][2])-int(mylist[0][1]))))
    else:
        length = [int(mylist[0][2])-int(mylist[0][1])]
        start = [0]
        for i in range(1,len(mylist)):
            length.append(int(mylist[i][2])-int(mylist[i][1]))
            start.append(int(mylist[i][1])-int(mylist[i-1][2]))
        lengthtmp = ','.join(map(str,length))+","
        starttmp = ','.join(map(str,start))+","
        out.write("{0}\t{1}\t{2}\t{3}\t0\t+\t{4}\t{5}\t0\t{6}\t{7}\t{8}\n".format(mylist[0][0],mylist[0][1],mylist[-1][2],barcode,mylist[0][1],mylist[-1][2],len(mylist),lengthtmp,starttmp))

def complex_cal(out):
    if len(out) >= 1000:
        return "F>1000"
    elif len(out) >= 500:
        return "F500_999"
    elif len(out) >= 200:
        return "F200_499"
    elif len(out) >= 100:
        return "F100_199"
    elif len(out) >= 50:
        return "F50_99"
    elif len(out) >= 40:
        return "F40_49"
    elif len(out) >= 30:
        return "F30_39"
    elif len(out) >= 20:
        return "F20_29"
    elif len(out) >= 15:
        return "F15_19"
    else:
        return "F{0}".format(len(out))

def deal(mylist,barcode,bed1,bed2,bed3):
    all_contain = []   # inter and intra
    IGVregion = []     # convert result to IGV bed
    intra_region = []  # just contain intra result
    complex_num = []   # BC complex number
    chrom = None       # chrom ID
    for line in mylist:
        if chrom != None and chrom != line[0]:
            tobed12(IGVregion,barcode,bed2)
            complex_num.append(complex_cal(intra_region))
            bed3.write(barcode+"\t"+str(len(intra_region))+"\t"+";".join(intra_region)+"\n")
            IGVregion = []
            intra_region = []
        chrom = line[0]
        all_contain.append("\t".join(line))
        IGVregion.append(line)
        intra_region.append("\t".join(line))
    tobed12(IGVregion,barcode,bed2)
    complex_num.append(complex_cal(intra_region))
    bed3.write(barcode+"\t"+str(len(intra_region))+"\t"+";".join(intra_region)+"\n")
    bed1.write(barcode+"\t"+str(len(all_contain))+"\t"+";".join(all_contain)+"\n")
    complex_num.append(complex_cal(all_contain))
    return complex_num

def fragment2BCbed(fin,outdir,prefix):
    all_complex_dict = {}
    intra_complex_dict = {}
    order = ('F1','F2','F3','F4','F5','F6','F7','F8','F9','F10','F11','F12','F13','F14','F15_19','F20_29','F30_39','F40_49','F50_99','F100_199','F200_499','F500_999','F>1000')
    for tmp_order in order:
        all_complex_dict[tmp_order] = 0
        intra_complex_dict[tmp_order] = 0
    outbed = os.path.join(outdir,prefix+".all.frag.bed")
    outbed2 = os.path.join(outdir,prefix+".IGV.bed")
    outbed3 = os.path.join(outdir,prefix+".samechrom.frag.bed")
    bed1 = open(outbed,'w')
    bed2 = open(outbed2,'w')
    bed3 = open(outbed3,'w')
    with open(fin,'r') as fl:
        region = []
        flag = None
        for line in fl:
            tmp = line.strip().split("-",1)
            barcode = tmp[0]
            if flag != None and flag != barcode:
                complex_num = deal(region,flag,bed1,bed2,bed3)
                all_complex_dict[complex_num[-1]] += 1
                for tmp_complex in complex_num[:len(complex_num)-1]:
                    intra_complex_dict[tmp_complex] += 1
                region = []
            flag = barcode
            region.append(tmp[1].split("\t")[0:3])
        complex_num = deal(region,flag,bed1,bed2,bed3)
        all_complex_dict[complex_num[-1]] += 1
        for tmp_complex in complex_num[:len(complex_num)-1]:
            intra_complex_dict[tmp_complex] += 1
    bed1.close()
    bed2.close()
    bed3.close()
    return(all_complex_dict,intra_complex_dict)


def convert(infile,allfile,PLECfile,PLISRSfile,fragment):
    PLEClist = []
    PLISRSlist = []
    for line in infile:
        tmp = line.strip().split("\t",2)
        if int(tmp[1]) < fragment:
            continue
        region = tmp[2].split(";")
        out = []
        for i in range(len(region)):
            for j in range(i,len(region)):
                score = 1/(2*int(tmp[1])-1)
                if i != j:
                    score *= 4
                pet1 = region[i].split()
                pet2 = region[j].split()
                pos1 = (int(pet1[1])+int(pet1[2]))//2
                pos2 = (int(pet2[1])+int(pet2[2]))//2
                PLEClist.append([0,pet1[0],pos1,0,0,pet2[0],pos2,1,score])
                out.append("\t".join(['0',pet1[0],str(pos1),'0','0',pet2[0],str(pos2),'1']))
        outrandom = set()
        for i in range(6):
            outrandom.add(random.choice(out))
        for i in outrandom:
            tmp = i.split("\t")
            PLISRSlist.append([tmp[0],tmp[1],int(tmp[2]),tmp[3],tmp[4],tmp[5],int(tmp[6]),tmp[7]])
    for i in sorted(PLEClist,key=itemgetter(1,5,2,6)):
        PLECfile.write("\t".join(map(str,i))+"\n")
        allfile.write("\t".join(map(str,i[:8]))+"\n")
    for i in sorted(PLISRSlist,key=itemgetter(1,5,2,6)):
        PLISRSfile.write("\t".join(map(str,i))+"\n")


def list2tree(mylist):
    resFrag={}
    for tmp in mylist:
        if tmp[0] in resFrag:
            tree = resFrag[tmp[0]]
            tree.add_interval(Interval(int(tmp[1]),int(tmp[2]), value={'name':tmp[3],'chr':tmp[0]}))
        else:
            tree = Intersecter()
            tree.add_interval(Interval(int(tmp[1]),int(tmp[2]), value={'name':tmp[3],'chr':tmp[0]}))
            resFrag[tmp[0]] = tree
    return resFrag

def overlap(resFrag,bed):
    myset = set()
    tmp = bed.split("\t")
    if tmp[0] in resFrag:
        result = resFrag[tmp[0]].find(int(tmp[1]),int(tmp[2]))
        for i in result:
            myset.add((i.value['chr'],i.start,i.end))
    # print(myset)
    return myset