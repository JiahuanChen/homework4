#! /usr/bin/python
"""
I worked with Haohan Gong.
Created on Mon Oct 24 15:05:43 2016
"""
import sys
from Bio import AlignIO
import math

#TASK 1: use Biopython to read file and save in alignment
alignment = AlignIO.read(sys.argv[1],"stockholm")

#TASK 2 calculate pi(x) and H(i)
#pi(x) in column_dist; H(i) in col_entropy

#pi(x)
seq_num = len(alignment)
seq_len = len(alignment[0])
column_dist = {}
#select a column
for col in range(0,seq_len):
    #initialize dictionary
    dist = {"A":0,"U":0,"G":0,"C":0,"-":0}
    #go throw the column and count
    for row in range(0,seq_num):
        dist[alignment[row][col]] += 1./seq_num
    column_dist[col] = dist

#H(i)
column_entropy = {}
for col in range(0,seq_len):
    h = 0;
    #h = -sigma(pi*log2(pi))
    for val in column_dist[col].values():
        if val != 0:
            h -= val*math.log(val,2)
    column_entropy[col] = h

#TASK3 pij(x,y)
nul = ["A","U","G","C","-"]
colpair_dist = {}
for i in range(0,seq_len-1):
    for j in range(i+1,seq_len):
        d = {}
        #initialize the dictionary, fill it with "AA":0, "AU":0 ...
        for nul1 in nul:
            for nul2 in nul:
                d[nul1+nul2] = 0

        for row in range(0,seq_num):
            d[alignment[row][i]+alignment[row][j]] += 1./seq_num
        colpair_dist[(i,j)] = d
                     
#TASK4 I(i,j)
colpair_mi = {}
#select i and j where i>j
for i in range(0,seq_len-1):
    for j in range(i+1,seq_len):
        mi = 0
	#sigma(sigma(pij*log2(pij/(pi*pj))))
        for nul1 in nul:
            for nul2 in nul:
                if colpair_dist[(i,j)][nul1+nul2] != 0:
                    pij = colpair_dist[(i,j)][nul1+nul2]
                    pi = column_dist[i][nul1]
                    pj = column_dist[j][nul2]
                    mi += pij*math.log(pij/(pi*pj),2)
        colpair_mi[(i,j)] = mi

#TASK5: print result
#item contance (key,value), where key is i = 0,1,2,...
#and value = entropy
item = column_entropy.items();
#change the position of keys and values
item = [[i[1],i[0]] for i in item];
#sort according the value in ascending order
item.sort()
for i in range(0,10):
    print(item[i][1])
    
# item contance (key,value), where key is (i,j) = (0,1),(0,2),...
#and value = mutual information, the rest part is like above
item = colpair_mi.items()
item = [[i[1],i[0]] for i in item];
#descending order
item.sort(reverse = True)
for i in range(0,50):
    print(str(item[i][1][0])+','+str(item[i][1][1]))
