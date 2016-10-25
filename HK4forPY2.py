#! /usr/bin/python
"""
Created on Mon Oct 24 15:05:43 2016

@author: blowv6
"""
from Bio import AlignIO
import math
#TASK 1 
#read file
alignment = AlignIO.read("tRNA.stock","stockholm")

#TASK 2
#pi(x) in column_dist; H(i) in col_entropy
seq_num = len(alignment)
seq_len = len(alignment[0])
column_dist = {}
for col in range(0,seq_len):
    dist = {"A":0,"U":0,"G":0,"C":0,"-":0}
    for row in range(0,seq_num):
        dist[alignment[row][col]] += 1./seq_num
    column_dist[col] = dist


column_entropy = {}
for col in range(0,seq_len):
    h = 0;
    for val in column_dist[col].values():
        if val != 0:
            #-sigma(pi*log2(pi))
            h -= val*math.log(val,2)
    column_entropy[col] = h

#TASK3 p i,j(x,y)
nul = ["A","U","G","C","-"]
colpair_dist = {}
for i in range(0,seq_len-1):
    for j in range(i+1,seq_len):
        d = {}
        #initialize the dictionary
        for nul1 in nul:
            for nul2 in nul:
                d[nul1+nul2] = 0
        for row in range(0,seq_num):
            d[alignment[row][i]+alignment[row][j]] += 1./seq_num
        colpair_dist[(i,j)] = d
                     
#TASK4 I(i,j)
colpair_mi = {}
for i in range(0,seq_len-1):
    for j in range(i+1,seq_len):
        mi = 0
        for nul1 in nul:
            for nul2 in nul:
                if colpair_dist[(i,j)][nul1+nul2] != 0:
                    pij = colpair_dist[(i,j)][nul1+nul2]
                    pi = column_dist[i][nul1]
                    pj = column_dist[j][nul2]
                    mi += pij*math.log(pij/(pi*pj),2)
        colpair_mi[(i,j)] = mi

#TASK5
item = column_entropy.items();
item = [[i[1],i[0]] for i in item];
item.sort()
for i in range(0,10):
    print(item[i][1])
    
item = colpair_mi.items()
item = [[i[1],i[0]] for i in item];
item.sort(reverse = True)
for i in range(0,50):
    print(item[i][1],item[i][0])