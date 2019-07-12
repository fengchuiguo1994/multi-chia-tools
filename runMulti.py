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

'''
function : deal with multi-chia data
version : 1.0.1
create : xingyu huang
begin : 2019-07-11
end : 
'''

parser = argparse.ArgumentParser(description = "deal with multi-chia data")
parser.add_argument("-i", "--input", required=True,help="the input bam/sam file")
parser.add_argument("-f", "--info", required=True,help="the input summary file")
parser.add_argument("-v", "--version", action = "version", version = "%(prog)s 1.0.1 by xyhuang")
parser.add_argument("-o", "--output", default = "output", help = "path of output file [output] ")
parser.add_argument("-p", "--prefix", default = "out", help = "prefix of output file [out] ")
parser.add_argument("-m", "--mapq", default = 30, type = int, help = "mapping quality [30] ")
parser.add_argument("-l", "--length", default = 50, type = int, help = "reference length [50] ")
parser.add_argument("-e", "--extend", default = 500, type = int, help = "extend length at 3' end [500] ")
parser.add_argument("-d", "--distance", default = 3000, type = int, help = "merge distance [3000] ")
parser.add_argument("-j", "--juicer", help = "path of Juicer tools")
parser.add_argument("-s", "--step", default = 1, type = int,help = "from step to end [1]")
parser.add_argument("-F", "--fragment", default = 2, type = int, help = "filter GEMs with the number of fragment < this value [2]")
parser.add_argument("-g", "--genome",help = "size of genome [hg19.chrom.size]")
parser.add_argument("-b", "--blacklist",help = "a file contain drop chromosome,one per line")
parser.add_argument("-a", "--anchor",help = "a bed file contain the chipseq anchor region")
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
        log.write("gems_detected : {0}\ntotal_reads : {1}\nmean_depth : {2}\npcr_duplication : {3}\n".format(infodict['gems_detected'],infodict['number_reads'],infodict['mean_depth'],infodict['pcr_duplication']))

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
        log.write("Rscript plot.r {0} {1}\n".format(args.output+"/"+args.prefix+".distance.txt",args.output+"/"+args.prefix+".distance.pdf"))
        print("Rscript {0}/plot.r {1} {2}".format(os.path.split(os.path.realpath(__file__))[0],args.output+"/"+args.prefix+".distance.txt",args.output+"/"+args.prefix+".distance.pdf"))
        returncode,returnresult = subprocess.getstatusoutput("Rscript plot.r {0} {1}".format(args.output+"/"+args.prefix+".distance.txt",args.output+"/"+args.prefix+".distance.pdf"))
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
            log.write("don't give juicer path and the genome size file,can't convert to hic file")
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
            log.write("don't give juicer path and the genome size file,can't convert to hic file")
        else:
            returncode,returnresult = subprocess.getstatusoutput("java -jar {0} pre -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000 {1} {2} {3} && java -jar {4} pre -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000 {5} {6} {7} && java -jar {8} pre -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000 {9} {10} {11}".format(args.juicer,allfile,allfile+".hic",args.genome,args.juicer,PLECfile,PLECfile+".hic",args.genome,args.juicer,PLISRSfile,PLISRSfile+".hic",args.genome))
            if returncode != 0:
                print ("[ERROR]: failed to generate HiC file : {0}\n".format(returnresult))
                exit()

        log.write("Step 5 use time : {0:.5f} s\n\n\n".format(time.time()-be))


    log.close()

'''
future work
## 1:use black list(drop some chromosome)
## 2:drop inter-interaction,just use intra-interaction(split complex to single and intra)
3:call cluster(use anchor)

plot work
'''