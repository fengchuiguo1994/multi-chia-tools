import networkx as nx
from networkx.algorithms import bipartite
from networkx.algorithms.bipartite import is_bipartite
from networkx.algorithms.bipartite import complete_bipartite_graph
from random import choice
import random
from bx.intervals.intersection import Intersecter, Interval
import sys
from collections import OrderedDict
import re

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
            myset.add(i.value['name'])
    # print(myset)
    return myset


anchor_file = '/public/home/xyhuang/xiaoqing/ChRDSeq/MH63_H3K4me3.anchor'
frag_file = '/public/home/xyhuang/multi-ChIA/MH63multi/combine/output/out.all.frag.bed'
cut_count = 3
prefix = "test"


anchorlist = []
anchor_map = OrderedDict()
with open(anchor_file,'r') as fin:
    for line in fin:
        tmp = line.strip().split()
        anchorlist.append(tmp)
        anchor_map[tmp[3]] = [tmp[:4]]
resFrag = list2tree(anchorlist)




countregion = {}
outregion = {}
with open(frag_file,'r') as fin:
    for line in fin:
        region_set = set()
        tmp = line.strip().split("\t",2)
        region = tmp[2].split(";")
        for i in region:
            tmpset = overlap(resFrag,i)
            region_set = region_set.union(tmpset)
        for j in region_set:
            if j not in countregion:
                countregion[j] = 0
            countregion[j] += 1
        if len(region_set) < 2: continue
        # aa = sorted(list(region_set),key=lambda x:(int(re.search("\d+",x[0])[0])))
        aa = sorted(list(region_set))
        for i in range(len(aa)-1):
            for j in range(i+1,len(aa)):
                if aa[i] not in outregion:
                    outregion[aa[i]] = {}
                if aa[j] not in outregion[aa[i]]:
                    outregion[aa[i]][aa[j]] = set()
                outregion[aa[i]][aa[j]].add(tmp[0])

fout_all = open(prefix+".colocalization.bed",'w')
fout = open(prefix+".co2.bed",'w')
for i in list(outregion.keys()):
    for j in list(outregion[i].keys()):
        if len(outregion[i][j]) < cut_count:
            del(outregion[i][j])
        else:
            fout.write("{0}\t{1}\t{2}\n".format(i,j,",".join(outregion[i][j])))
            fout_all.write("{0}\t{1}\t{2}\n".format(i,j,",".join(outregion[i][j])))
        if len(outregion[i]) == 0:
            del(outregion[i])
co_region = outregion.copy()
fout.close()


N = 3
while len(outregion) > 0:
    fout = open("{0}.co{1}.bed".format(prefix,N),'w')
    for i in list(outregion.keys()):
        for j in list(outregion[i].keys()):
            if len(outregion[i][j]) < cut_count:
                del(outregion[i][j])
            if len(outregion[i]) == 0:
                del(outregion[i])
            
            if j in co_region:
                for k in list(co_region[j].keys()):
                    outregion[i+"\t"+j] = {}
                    outregion[i+"\t"+j][k] = outregion[i][j].intersection(co_region[j][k])
                    
    for i in list(outregion.keys()):
        if len(i.split("\t")) < N-1:
            del(outregion[i])
            continue
        elif len(i.split("\t")) > N-1:
            print("some problem")
            exit()
        for j in list(outregion[i].keys()):
            if len(outregion[i][j]) < cut_count:
                del(outregion[i][j])
            else:
                fout.write("{0}\t{1}\t{2}\n".format(i,j,",".join(outregion[i][j])))
                fout_all.write("{0}\t{1}\t{2}\n".format(i,j,",".join(outregion[i][j])))
            if len(outregion[i]) == 0:
                del(outregion[i])
    fout.close()
    N += 1
fout_all.close()