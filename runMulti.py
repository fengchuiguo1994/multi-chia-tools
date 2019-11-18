import argparse
import time
import subprocess
import sys
import pysam
import re
import os
import pybedtools
import random
import function
from scipy.stats import hypergeom
import pandas as pd
from numpy import log10

'''
function : deal with multi-chia data
version : 2.0.0
create : xingyu huang
begin : 2019-10-11
end : 
'''

parser = argparse.ArgumentParser(description = "deal with multi-chia data")
parser.add_argument("-i", "--input", required=True,help="the input bam/sam file")
parser.add_argument("-f", "--info", required=True,help="the input summary file")
parser.add_argument("-v", "--version", action = "version", version = "%(prog)s 2.0.0 by xyhuang")
parser.add_argument("-j", "--juicer", help = "the path of Juicer tools,if set,the -g is must")
parser.add_argument("-g", "--genome",help = "size of genome file;like hg19.chrom.size")
parser.add_argument("-o", "--output", default = "output", help = "path of output file [default:output] ")
parser.add_argument("-p", "--prefix", default = "out", help = "prefix of output file [default:out] ")
parser.add_argument("-m", "--mapq", default = 30, type = int, help = "mapping quality [default:30] ")
parser.add_argument("-l", "--length", default = 50, type = int, help = "the length of reads match to reference genome [default:50] ")
parser.add_argument("-e", "--extend", default = 500, type = int, help = "extend length at 3' end(from 5' to 3') [default:500] ")
parser.add_argument("-d", "--distance", default = 3000, type = int, help = "merge two fragment if the distance is less 3000 [default:3000] ")
parser.add_argument("-s", "--step", default = 1, type = int,help = "from the step to end [default:1]")
parser.add_argument("-F", "--fragment", default = 2, type = int, help = "filter GEMs with the number of fragment < this value [default:2]")
parser.add_argument("-b", "--blacklist",help = "a file contain drop chromosome,one chromosome per line")
parser.add_argument("-a", "--anchor",help = "a bed file contain the chipseq anchor region,we will use it to find cluster")
args = parser.parse_args()



if __name__ == "__main__":
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    logfile = os.path.join(args.output,args.prefix+"_log.txt")
    if args.step == 1:
        log = open(logfile,'w')
    else:
        log = open(logfile,'a')


#  get summary file infomation   #
    if args.step == 1:
        args.step += 1
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + " Step 1: get align static ...")
        log.write(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + " Step 1: get align static ...\n")
        be = time.time()

        fin = open(args.info,'r').readlines()
        infodict = {}
        flag = fin[0].split(",")
        val = fin[1].split(",")
        for i,j in zip(flag,val):
            infodict[i] = j
        log.write("gems_detected : {0}\nmean_dna_per_gem : {1}\nmolecule_length_mean : {2}\nmolecule_length_stddev : {3}\ntotal_reads : {4}\nmedian_insert_size : {5}\nmean_depth : {6}\nzero_coverage : {7}\npcr_duplication : {8}\n".format(infodict['gems_detected'],infodict['mean_dna_per_gem'],infodict['molecule_length_mean'],infodict['molecule_length_stddev'],infodict['number_reads'],infodict['median_insert_size'],infodict['mean_depth'],infodict['zero_coverage'],infodict['pcr_duplication']))

        log.write("Step 1 use time : {0:.5f} s\n\n\n".format(time.time()-be))

#  convert bam file to bed file   #
    if args.step == 2:
        args.step += 1
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + " Step 2: convert bam file to bed file ...")
        log.write(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + " Step 2: convert bam file to bed file ...\n")
        be = time.time()

        black = set()
        if args.blacklist == None:
            log.write("we will not drop any chromosome\n")
        else:
            with open(args.blacklist,'r') as fin:
                for line in fin:
                    line = line.strip()
                    black.add(line)
            log.write("we will drop these chromosome : {0}\n".format(",".join(black)))

        samfile = function.readfile(args.input)
        total_reads,nn,typeBC,R1total,R2total,R1barcode,R2barcode,R1countmap,R2countmap,R1uniqmap,R2uniqmap,R1Q30map,R2Q30map,R1lenmap,R2lenmap,R1black,R2black = function.bam2bed(samfile,args.output,args.prefix,args.mapq,args.length,args.extend,black)
        with open(os.path.join(args.output,args.prefix+".BC.txt"),'w') as fout:
            for tmpkey in typeBC.keys():
                fout.write(tmpkey+"\t"+str(typeBC[tmpkey])+"\n")
        del(typeBC)
        log.write("reads_number : {0}\ntotal_reads : {1}\nR1 total_reads : {2}\nR2 total_reads : {3}\nR1 with barcode : {4}\nR2 with barcode : {5}\nR1 with barcode,map : {6}\nR2 with barcode,map : {7}\nR1 with barcode,map,uniq : {8}\nR2 with barcode,map,uniq : {9}\nR1 with barcode,map,uniq,Q{10} : {11}\nR2 with barcode,map,uniq,Q{12} : {13}\nR1 with barcode,map,uniq,Q{14},len{15} : {16}\nR2 with barcode,map,uniq,Q{17},len{18} : {19}\nR1 with barcode,map,uniq,Q{20},len{21},chrom : {22}\nR2 with barcode,map,uniq,Q{23},len{24},chrom : {25}\n".format(total_reads,nn,R1total,R2total,R1barcode,R2barcode,R1countmap,R2countmap,R1uniqmap,R2uniqmap,args.mapq,R1Q30map,args.mapq,R2Q30map,args.mapq,args.length,R1lenmap,args.mapq,args.length,R2lenmap,args.mapq,args.length,R1black,args.mapq,args.length,R2black))

        log.write("Step 2 use time : {0:.5f} s\n\n\n".format(time.time()-be))

#  bed file sort and merge        #
    if args.step == 3:
        args.step += 1
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + " Step 3: bed file sort and merge ...")
        log.write(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + " Step 3: bed file sort and merge ...\n")
        be = time.time()

        bedfile = pybedtools.BedTool(os.path.join(args.output,args.prefix+".cv.bed"))
        outbedfile = os.path.join(args.output,args.prefix+".merge.bed")
        merge = bedfile.sort().merge(d = args.distance,c=4,o="count,collapse")
        merge.saveas(outbedfile)
        log.write("bed_record : {0}\nfragment_record : {1}\n".format(str(bedfile.count()),str(merge.count())))

        log.write("Step 3 use time : {0:.5f} s\n\n\n".format(time.time()-be))
    

#  convert fragment to BCbed    #
    if args.step == 4:
        args.step += 1
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + " Step 4: convert fragment to BCbed ...")
        log.write(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + " Step 4: convert fragment to BCbed ...\n")
        be = time.time()

        bedfile = os.path.join(args.output,args.prefix+".merge.bed")
        all_complex_dict,intra_complex_dict = function.fragment2BCbed(bedfile,args.output,args.prefix)
        log.write("for all connection\n")
        out = []
        order = ('F1','F2','F3','F4','F5','F6','F7','F8','F9','F10','F11','F12','F13','F14','F15_19','F20_29','F30_39','F40_49','F50_99','F100_199','F200_499','F500_999','F>1000')
        total_complex1 = total_complex2 = 0
        for i in order:
            out.append("{0} : {1}".format(i,all_complex_dict[i]))
            total_complex1 += all_complex_dict[i]
            if i != "F1":
                total_complex2 += all_complex_dict[i]

        log.write("\n".join(out)+"\n")
        log.write("chromatin complexes: {0}\nchromatin complexes >=2 : {1}\n\n\n".format(total_complex1,total_complex2))
        
        log.write("for just same chromosome connection\n")
        out = []
        total_complex1 = total_complex2 = 0
        for i in order:
            out.append("{0} : {1}".format(i,intra_complex_dict[i]))
            total_complex1 += intra_complex_dict[i]
            if i != "F1":
                total_complex2 += intra_complex_dict[i]

        log.write("\n".join(out)+"\n")
        log.write("chromatin complexes: {0}\nchromatin complexes >=2 : {1}\n".format(total_complex1,total_complex2))
        
        with open(args.output+"/"+args.prefix+".IGV.bed",'r') as fin,open(args.output+"/"+args.prefix+".distance.txt",'w') as fout:
            for line in fin:
                tmp = line.strip().split()
                if int(tmp[9]) > 1:
                    for i in tmp[-1].split(","):
                        if i != "0" and i != "":
                            fout.write("{0}\t{1}\n".format(i,tmp[9]))
        log.write("Rscript {0}/R/distance.r {1} {2}\n".format(os.path.split(os.path.realpath(__file__))[0],args.output+"/"+args.prefix+".distance.txt",args.output+"/"+args.prefix+".distance.pdf"))
        print("Rscript {0}/R/distance.r {1} {2}\n".format(os.path.split(os.path.realpath(__file__))[0],args.output+"/"+args.prefix+".distance.txt",args.output+"/"+args.prefix+".distance.pdf"))
        returncode,returnresult = subprocess.getstatusoutput("Rscript {0}/R/distance.r {1} {2}".format(os.path.split(os.path.realpath(__file__))[0],args.output+"/"+args.prefix+".distance.txt",args.output+"/"+args.prefix+".distance.pdf"))
        if returncode != 0:
            print ("[ERROR]: failed to plot : {0}\n".format(returnresult))
            exit()

        log.write("Step 4 use time : {0:.5f} s\n\n\n".format(time.time()-be))


#  convert fragment to HiC    #
    if args.step == 5:
        args.step += 1
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + " Step 5: convert fragment to HiC ...")
        log.write(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + " Step 5: convert fragment to HiC ...\n")
        be = time.time()

        infile = os.path.join(args.output,args.prefix+".all.frag.bed")
        allfile = os.path.join(args.output,args.prefix+".total.all.bedpe")
        PLECfile = os.path.join(args.output,args.prefix+".total.PLEC.bedpe")
        PLISRSfile = os.path.join(args.output,args.prefix+".total.PLISRS.bedpe")
        infl = open(infile,'r')
        allout = open(allfile,'w')
        PLECout = open(PLECfile,'w')
        PLISRSout = open(PLISRSfile,'w')
        function.convert(infl,allout,PLECout,PLISRSout,args.fragment)
        infl.close()
        allout.close()
        PLECout.close()
        PLISRSout.close()
        if args.juicer == None or args.genome == None:
            log.write("don't give juicer path and the genome size file,can't convert to hic file\n")
        else:
            returncode,returnresult = subprocess.getstatusoutput("java -jar {0} pre -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000 {1} {2} {3} && java -jar {4} pre -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000 {5} {6} {7} && java -jar {8} pre -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000 {9} {10} {11}".format(args.juicer,allfile,allfile+".hic",args.genome,args.juicer,PLECfile,PLECfile+".hic",args.genome,args.juicer,PLISRSfile,PLISRSfile+".hic",args.genome))
            if returncode != 0:
                print ("[ERROR]: failed to generate HiC file : {0}\n".format(returnresult))
                exit()
            
        infile = os.path.join(args.output,args.prefix+".samechrom.frag.bed")
        allfile = os.path.join(args.output,args.prefix+".intra.all.bedpe")
        PLECfile = os.path.join(args.output,args.prefix+".intra.PLEC.bedpe")
        PLISRSfile = os.path.join(args.output,args.prefix+".intra.PLISRS.bedpe")
        infl = open(infile,'r')
        allout = open(allfile,'w')
        PLECout = open(PLECfile,'w')
        PLISRSout = open(PLISRSfile,'w')
        function.convert(infl,allout,PLECout,PLISRSout,args.fragment)
        infl.close()
        allout.close()
        PLECout.close()
        PLISRSout.close()
        if args.juicer == None or args.genome == None:
            log.write("don't give juicer path and the genome size file,can't convert to hic file\n")
        else:
            returncode,returnresult = subprocess.getstatusoutput("java -jar {0} pre -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000 {1} {2} {3} && java -jar {4} pre -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000 {5} {6} {7} && java -jar {8} pre -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000 {9} {10} {11}".format(args.juicer,allfile,allfile+".hic",args.genome,args.juicer,PLECfile,PLECfile+".hic",args.genome,args.juicer,PLISRSfile,PLISRSfile+".hic",args.genome))
            if returncode != 0:
                print ("[ERROR]: failed to generate HiC file : {0}\n".format(returnresult))
                exit()

        log.write("Step 5 use time : {0:.5f} s\n\n\n".format(time.time()-be))


#  convert fragment to ChIA-PET    #
    if args.step == 6:
        args.step += 1
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + " Step 6: convert fragment to ChIA-PET ...")
        log.write(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + " Step 6: convert fragment to ChIA-PET ...\n")
        be = time.time()

        if args.anchor == None:
            log.write("there no anchor file,so we will use the singleton data to generate the high coverage region instead anchor\n")
            '''
            we will do the work in next version
            '''
        else:
            N = 0 # fragment anchor count
            Ntmp = 0 # fragment count
            log.write("we will build a anchor region tree\n")
            anchorlist = []
            with open(args.anchor,'r') as fin:
                for line in fin:
                    tmp = line.strip().split()
                    anchorlist.append(tmp)
            resFrag = function.list2tree(anchorlist)
            # for i in resFrag.keys():
                # print(i)
            
            infile = os.path.join(args.output,args.prefix+".all.frag.bed")
            countregion = {}
            # outregion = []
            outregion = {}
            with open(infile,'r') as fin:
                for line in fin:
                    region_set = set()
                    tmp = line.strip().split("\t",2)
                    region = tmp[2].split(";")
                    for i in region:
                        N += 1
                        Ntmp += 1
                        tmpset = function.overlap(resFrag,i)
                        if len(tmpset) > 1:
                            N += len(tmpset)-1
                        region_set = region_set.union(tmpset)
                    for i in region_set:
                        region = i[0]+"\t"+str(i[1])+"\t"+str(i[2])
                        if region not in countregion:
                            countregion[region] = 0
                        countregion[region] += 1
                    if len(region_set) < 2: continue
                    # aa = list(region_set)
                    aa = sorted(list(region_set),key=lambda x:(int(re.search("\d+",x[0])[0]),x[1]))
                    for i in range(len(aa)-1):
                        for j in range(i+1,len(aa)):
                            # print(aa[i][0]+"\t"+str(aa[i][1])+"\t"+str(aa[i][2])+"\t"+aa[j][0]+"\t"+str(aa[j][1])+"\t"+str(aa[j][2]))
                            # outregion.append([aa[i][0]+"\t"+str(aa[i][1])+"\t"+str(aa[i][2]),aa[j][0]+"\t"+str(aa[j][1])+"\t"+str(aa[j][2])])
                            tmpregion = (aa[i][0]+"\t"+str(aa[i][1])+"\t"+str(aa[i][2]),aa[j][0]+"\t"+str(aa[j][1])+"\t"+str(aa[j][2]))
                            if tmpregion not in outregion:
                                outregion[tmpregion] = 0
                            outregion[tmpregion] += 1
            chiafile = os.path.join(args.output,args.prefix+".chiapet.cluster")
            chiaout = open(chiafile,'w')
            for i,j in outregion.keys():
                chiaout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i,j,outregion[(i,j)],countregion[i],countregion[j],N))
            chiaout.close()
        log.write("fragment count:{0}\nfragment anchor count:{1}\n".format(Ntmp,N))

        # use hypergeom to calculator Significant cluster
        chiapvalue = os.path.join(args.output,args.prefix+".chiapet.cluster.withpvalue")
        chiaFDR = os.path.join(args.output,args.prefix+".chiapet.cluster.withpvalue.FDR.txt")
        chialist = []
        chiapd = pd.read_table(chiafile,header=None)
        chiapd.columns=["chr1",'start1','end1',"chr2",'start2','end2','petcount','coverage1','coverage2','totalaln']
        chiapd['pvalue'] = hypergeom.sf(chiapd['petcount'],chiapd['totalaln'],chiapd['coverage1'],chiapd['coverage2'])

        # calculator FDR
        from rpy2.robjects.packages import importr
        from rpy2.robjects.vectors import FloatVector
        from rpy2.robjects import pandas2ri
        stats = importr('stats')
        
        pVal = chiapd['pvalue']
        fdr = stats.p_adjust(FloatVector(pVal),method='BH')
        fdrPD = pandas2ri.ri2py_vector(fdr)
        # chiapd.insert(10,'FDR',fdrPD)
        chiapd['FDR'] = fdrPD

        # pv = chiapd['pvalue']
        # pvlist = [i for i in pv]
        # lengthpv = len(pvlist)
        # pvsort = [i for i in pvlist]
        # pvsort.sort()
        # pvFDR = []
        # for i in pvlist:
        #     qv = i*lengthpv/(pvsort.index(i)+1)
        # pvFDR.append(qv)
        # chiapd['FDR']= pvFDR

        chiapd['-log10pvalue'] = -log10(chiapd['pvalue'])
        chiapd['-log10qvalue'] = -log10(chiapd['FDR'])

        chiapd['-log10pvalue'] = chiapd['-log10pvalue'].map('{:.4f}'.format)
        chiapd['-log10qvalue'] = chiapd['-log10qvalue'].map('{:.4f}'.format)
        chiapd['pvalue'] = chiapd['pvalue'].map('{:.4e}'.format)
        chiapd['FDR'] = chiapd['FDR'].map('{:.4e}'.format)
        chiapd.to_csv(chiapvalue,sep="\t",index=False)

        with open(chiapvalue,'r') as fin,open(chiaFDR,'w') as fout:
            for line in fin:
                if line.strip().endswith("qvalue"):
                    fout.write(line)
                else:
                    tmp = line.strip().split()
                    if int(tmp[6]) > 1 and float(tmp[11]) < 0.01:
                        fout.write(line)

        log.write("Step 6 use time : {0:.5f} s\n\n\n".format(time.time()-be))

    log.close()

'''
future work
## 1:use black list(drop some chromosome)
## 2:drop inter-interaction,just use intra-interaction(split complex to single and intra)
## 3:call cluster(use anchor)

plot work
'''
'''
from bx.intervals.intersection import Intersecter, Interval

resFrag={}
with open("/public/home/xyhuang/xiaoqing/ChRDSeq/MH63_H3K4me3.anchor",'r') as fin:
    for line in fin:
        tmp = line.strip().split("\t")
        if tmp[0] in resFrag:
            tree = resFrag[tmp[0]]
            tree.add_interval(Interval(int(tmp[1]),int(tmp[2]), value={'name':tmp[3],'chr':tmp[0]}))
        else:
            tree = Intersecter()
            tree.add_interval(Interval(int(tmp[1]),int(tmp[2]), value={'name':tmp[3],'chr':tmp[0]}))
            resFrag[tmp[0]] = tree

def overlap(resFrag,bed):
    myset = set()
    tmp = bed.split("\t")
    if tmp[0] in resFrag:
        result = resFrag[tmp[0]].find(int(tmp[1]),int(tmp[2]))
        for i in result:
            myset.add((i.value['chr'],i.start,i.end))
    return myset

with open("/public/home/xyhuang/multi-ChIA/MH63multi/RMCD001/RMCD001/outs/output/out.all.frag.bed",'r') as fin:
    for line in fin:
        region_set = set()
        tmp = line.strip().split("\t",2)
#        print(tmp[2])
        region = tmp[2].split(";")
        for i in region:
            region_set = region_set.union(overlap(resFrag,i))
        print(region_set)
'''