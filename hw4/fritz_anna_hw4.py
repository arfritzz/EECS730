#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import csv

# 1. read string from file 
# seq is the input string 
f1 = open('seq.fa', 'r')
seq = ""
for line in f1:
    if line[0] != '>': 
        seq += line
f1.close()
seq = seq.replace("\n","")
#print(seq)

# 2. read in table from file 
always_print = False
counter = 0
mod_file = open("hmm_trimmed.txt", "w") 

f1 = open('AAA.hmm', 'r')
for line in f1:
    if '//' in line: 
        always_print = False
    if always_print or 'HMM  ' in line: 
        counter+=1
        always_print = True
        mod_file.write(line)
        # print(line.strip().split("\t"))

f1.close()
mod_file.close()

with open("hmm_trimmed.txt", 'r') as f:
    LoL=[x.split() for x in f]
        
# data is the hmm data 
data = pd.DataFrame(LoL)

# this take is the state to letter conversion 
state_to_letter = data.iloc[: , 21:]
state_to_letter = state_to_letter.replace(to_replace='None', value=np.nan).dropna()
state_to_letter = state_to_letter.reset_index(drop=True)
state_to_letter.columns=['A','B','C','D','E']

#delete the first row and make it the index 
data.columns = data.iloc[0]
data = data.drop(data.index[0])

#this is the table with the emission frequencies for the match
e_M = data[1::3]
e_M = e_M.dropna(axis=1)
e_M = e_M.iloc[: , 1:]
e_M = e_M.reset_index(drop=True)
e_M = e_M.astype(float)

#this is the table with the emission frequencies for the insert
e_I = data[2::3]
e_I = e_I.shift(periods=1, axis="columns")
e_I = e_I.dropna(axis=1)
e_I = e_I.reset_index(drop=True)
e_I = e_I.astype(float)

# this is the table with the transition frequencies 
Trans = data[::3]
Trans = Trans.dropna(axis=1)
Trans.columns=Trans.iloc[0]
Trans = Trans.drop(Trans.index[0])
Trans = Trans.reset_index(drop=True)
Trans = Trans.replace('*','0')
Trans = Trans.astype(float)





# # helpful links 
# 
# http://www2.cs.uh.edu/~ceick/ML/HMM_in_BI.pdf
# 
# https://www.cis.upenn.edu/~cis262/notes/Example-Viterbi-DNA.pdf

# # steps 
# 
# 1. make a table for Vm, Vi, Vd
#     - fill match table, then insert table, then delete table
# 2. traceback table 

#these are the v matrix 
v_M = pd.DataFrame(np.zeros((133, len(seq)+1)))
v_D = pd.DataFrame(np.zeros((133, len(seq)+1)))
v_I = pd.DataFrame(np.zeros((133, len(seq)+1)))

# these hold the directions for where each value came from 
d_M = pd.DataFrame(np.zeros((133, len(seq)+1)))
d_D = pd.DataFrame(np.zeros((133, len(seq)+1)))
d_I = pd.DataFrame(np.zeros((133, len(seq)+1)))



# initalizing 
v_M[0][0] = 1
v_D[0][0] = 1
v_I[0][0] = 1

# filling the tables 

for i in range(1,len(v_M)): 
    for j in range(1,len(seq)+1): 
        current_letter = seq[j-1]
        
        
        #MATCH MATRIX 
        v_max = max((v_M.loc[i-1,j-1] - Trans['m->m'][i]), (v_I.loc[i-1,j-1] - Trans['i->m'][i]), (v_D.loc[i-1,j-1] - Trans['d->m'][i]))
            
        v_M.loc[i,j] = -e_M[current_letter][i] + e_M[current_letter][0] + v_max
        
        # direction for value from M matrix
        if (v_max == (v_M.loc[i-1,j-1] - Trans['m->m'][i])): 
            d_M.loc[i,j] = 'M'
        elif (v_max == (v_I.loc[i-1,j-1] - Trans['i->m'][i])):
            d_M.loc[i,j] = 'I'
        else: 
            d_M.loc[i,j] = 'D'
            
        
        #DELETE MATRIX
        v_max =  max((v_M.loc[i-1,j] - Trans['m->d'][i]), (v_D.loc[i-1,j] - Trans['d->d'][i]))
        
        # direction for value from D matrix 
        if (v_max == (v_M.loc[i-1,j] - Trans['m->d'][i])):
            d_D.loc[i,j] = 'M'
        else: 
            d_D.loc[i,j] = 'D'
        
        v_D.loc[i,j] = v_max
            
            
            
        #INSERT MATRIX 
        v_max =  max((v_M.loc[i,j-1] - Trans['m->i'][i]), (v_I.loc[i,j-1] - Trans['i->i'][i]))
        
        # direction for value from I matrix
        if (v_max == (v_M.loc[i,j-1] - Trans['m->i'][i])):
            d_I.loc[i,j] = 'M'
        else: 
            d_I.loc[i,j] = 'I'
        
        v_I.loc[i,j] = -e_I[current_letter][i] + e_I[current_letter][0] + v_max


#d_M.to_csv('dm.csv')
#v_M.to_csv('vm.csv')

#these will be the final sequences 
seq1 = []
seq_mid = []
seq2 = []

# get max value and print as score 

i = 132
#i = 133
j = len(seq)

current = max(v_M.loc[i,j], v_D.loc[i,j], v_I.loc[i,j])

if (current == v_M.loc[i,j]): 
    direction = d_M.loc[i,j]
elif (current == v_D.loc[i,j]): 
    direction = d_D.loc[i,j]
else: 
    direction = d_I.loc[i,j]


# PRINTING
# 1 prints score 
# 2 looks at the table and prints the state that corresponds to the current letter 
#   then prints the sequences 
print('score: ' + str(current))

while(i>0 and j>0):

    current_letter = seq[j-1]
     
    #what to print 
    if (direction == 'M'):
        seq1 += state_to_letter['B'][i-1]
        if (-e_M[current_letter][i] + e_M[current_letter][0] > 0):
            seq_mid += '+'
        else: 
            seq_mid += ' '
        seq2 += seq[j-1] 
        i = i-1
        j = j-1
        
        direction = d_M.loc[i,j]
    elif (direction == 'D'):
        seq1 += state_to_letter['B'][i-1]
        seq_mid += ' '
        seq2 += '-'
        i = i-1
        direction = d_D.loc[i,j]

    else: 
        seq1 += state_to_letter['B'][i-1]
        if (-e_I[current_letter][i] + e_I[current_letter][0] > 0):
            seq_mid += '+'
        else: 
            seq_mid += ' '
        seq2 += seq[j-1]
        j=j-1
        direction = d_I.loc[i,j]

# this prints if the front of the sequence has been deleted 
# prof said not to do this 
#while (j>0): 
#    seq2 += seq[j-1]
#    seq1 += '-'
#    seq_mid += ' '
#    j = j -1

seq1.reverse()
seq2.reverse()
seq_mid.reverse()

seq1ToStr = ' '.join([str(elem) for elem in seq1])
seq2ToStr = ' '.join([str(elem) for elem in seq2])
seqMToStr = ' '.join([str(elem) for elem in seq_mid])

print(seq1ToStr)
print(seqMToStr)
print(seq2ToStr)
