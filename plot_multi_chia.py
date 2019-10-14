# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 09:17:27 2019

@author: huang
"""

import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle

# import matplotlib.ticker as ticker
# import matplotlib.patches as patches
# import numpy as np

import pandas as pd
import re

def getmax(df):
    return max(df[1],df[2],df[4],df[5])
def getmin(df):
    return min(df[1],df[2],df[4],df[5])

def readfile(infile):
    data=pd.read_table(infile,sep="\t",header=None)
    if len(data.columns) == 6:
        data["high"] = 2
    else:
        data["high"] = data[6]
    data["start"] = data.apply(getmin,axis = 1)
    data["end"] = data.apply(getmax,axis = 1)
    return data

def readmulti(multi):
    multiregion = []
    with open(multi,'r') as fin:
        for line in fin:
            tmp = line.strip().split("\t",maxsplit=2)
            for i in tmp[2].split(";"):
                ll = i.split("\t")
                multiregion.append([tmp[0],int(tmp[1]),ll[0],int(ll[1]),int(ll[2])])
    return(pd.DataFrame(multiregion))

def region_split(string):
    chrom,start,end = re.split(r":|-",region)
    return chrom,int(start),int(end)

def plot_bedpe(df,axs):
    before = int((df[1]+df[2])/2)
    after = int((df[4]+df[5])/2)
    min2 = int((before+after)/2)
    pp1 = mpatches.PathPatch(
    Path([(before, 0), (min2,2*df['high']), (after, 0),(before, 0)],
         [Path.MOVETO, Path.CURVE3, Path.CURVE3,Path.CLOSEPOLY]),
    fc="none", transform=axs.transData,color="red")
    axs.add_patch(pp1)

def plot_multi(df,axs,drop=True):
    flag = None
    out = []
    tmp = []
    for index,row in df.iterrows():
        if flag != None and flag != row['barcode']:
            out.append(tmp)
            tmp = []
        flag = row['barcode']
        tmp.append((int(row['start']),int(row['end'])))
    out.append(tmp)
    outregion = []
    if drop == True:
        for line in out:
            if len(line) == 1:
                continue
            tt = min(line, key=lambda x: x[1])[1]
            line.append(max(line, key=lambda x: x[0])[0])
            line.append(tt)
            outregion.append(list(reversed(line)))
    else:
        for line in out:
            tt = min(line, key=lambda x: x[1])[1]
            line.append(max(line, key=lambda x: x[0])[0])
            line.append(tt)
            outregion.append(list(reversed(line)))
    outregion.sort(key=lambda x: (x[0],x[1]))
    # print(outregion)
    
    flagregion = []
    for oo,i in enumerate(outregion):
        flag_overlap = False
        if len(flagregion) == 0:
            flagregion.append([i[0],i[1]])
            for j in i[2:]:
                rect = Rectangle((j[0],1-0.2),j[1]-j[0],0.4, color="red", fill=True)
                axs.add_patch(rect)
            axs.plot([i[0],i[1]],[1,1],color='black', linestyle='-')
        else:
            for n,m in enumerate(flagregion):
                if (((i[0] > m[0]) & (i[0] < m[1])) | ((i[1] > m[0]) & (i[1] < m[1])) | ((i[0] < m[0]) & (i[1] > m[1]))):
                    pass
                else:
                    flag_overlap = n + 1
                    flagregion[n][1] = i[1]
                    break
            if flag_overlap == False:
                flagregion.append([i[0],i[1]])
                flag_overlap = len(flagregion)
            for j in i[2:]:
                rect = Rectangle((j[0],flag_overlap-0.2),j[1]-j[0],0.4, color="red", fill=True)
                axs.add_patch(rect)
            axs.plot([i[0],i[1]],[flag_overlap,flag_overlap],color='black', linestyle='-')
    return len(flagregion)+1

if __name__ == '__main__':
    # get interaction file
    infile = r"C:\Users\huang\Desktop\cluters.txt"
    bedpe = readfile(infile)
    bedpe_group = bedpe.groupby(0)
    # bedpe_group.get_group('chr1')
    # get multi-chia file
    infile = r"C:\Users\huang\Desktop\out.samechrom.frag.bed"
    multi = readmulti(infile)
    multi.columns = ["barcode","count","chrom","start","end"]
    multi_group = multi.groupby("chrom")
    # multi_group.get_group('chr1')
    
    # give a plot region and filter data
    region = "chr1:1500000-1600000"
    chrom,start,end = region_split(region)
    tmp = bedpe_group.get_group(chrom)
    flt_bedpe = tmp[(tmp['start'] > start) & (tmp['end'] < end)]
    tmp = multi_group.get_group(chrom)
    # flt_multi = tmp[((tmp['start'] > start) & (tmp['start'] < end)) | ((tmp['end'] > start) & (tmp['end'] < end)) | ((tmp['start'] < start) & (tmp['end'] > end))].sort_values(by = ["start","end"],ascending = [True,True])
    flt_multi = tmp[((tmp['start'] > start) & (tmp['start'] < end)) | ((tmp['end'] > start) & (tmp['end'] < end)) | ((tmp['start'] < start) & (tmp['end'] > end))]
    del(tmp)
    
    # plot bedpe
    Path = mpath.Path
    # fig = plt.figure(figsize=(4,8))
    # axs = fig.add_subplot(211)
    axs = plt.subplot(211)
    flt_bedpe.apply(plot_bedpe,axis = 1,axs=axs)
    axs.set_title('chia-pet cluster')
    # axs.set_xlabel(chrom)
    axs.set_ylabel('petcount')
    axs.axis([start,end,2,8])
    axs.set_xticklabels([])
    
    # axs1 = fig.add_subplot(212)
    axs1 = plt.subplot(212)
    aa = plot_multi(flt_multi,axs=axs1,drop=True) # 去掉singleton的。
    axs1.set_title('multi-chia reads linker')
    axs1.set_xlabel(chrom)
    axs1.set_ylabel('barcode')
    axs1.axis([start,end,0,aa])

    # plt.show()
    plt.savefig("D:/temp.pdf",dpi=300)