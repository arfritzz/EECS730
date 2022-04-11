#!/usr/bin/python

#from Bio import SeqIO
import sys
import array as arr

f1=open("seq2.fa", "r")
# this will be the first line which is the name of the sequence 
f1.readline()
seq1 = list(f1.readline())
#print(seq1)
f1.close()

f1=open("seq1.fa", "r")
# this will be the first line which is the name of the sequence 
f1.readline()
seq2 = list(f1.readline())
#print(seq1)
f1.close()

#print(seq1, seq2)

#created matrix a and initalized all values to 0 
a = [[0 for x in range(len(seq2)+1)] for x in range(len(seq1)+1)]
d = [["A" for x in range(len(seq2)+1)] for x in range(len(seq1)+1)]

#for i in range(len(a)):
#    for j in range(len(a[i])):
#        print((i,j, a[i][j]) , end=" ")
#    print()

# initialize parameters for best location to
# know where to start when reading back 
best = 0
best_cord = (0,0)

gap = -3


#fn to match if the sequences are the same  
def match (i,j):
    #print((seq1[i-1],seq2[j-1]))
    if (seq1[i-1] == seq2[j-1]): 
        return(1)
    else: 
        return(-2)


# for loop to fill values in a based on alignment algorithm 
for i in range(1,len(a)):
    for j in range(1,len(a[i])):
        #print((i,j,a[i][j]),end=' ')
        a[i][j] = max(a[i][j-1]+ gap, 
                      a[i-1][j] + gap,
                      a[i-1][j-1] + match(i,j))
        if ( seq1[i-1] == seq2[j-1]): 
            # came from left
            if (a[i][j] == a[i][j-1]+ gap):
                d[i][j] = "L"
            #diagonal
            elif (a[i][j] == a[i-1][j-1]+ match(i,j)):
                d[i][j] = "D"
            # came from up 
            else:
                d[i][j] = "U"
        # always start with a diagonal? 
        elif (i == len(a)-1 and j == len(a[i])-1): 
            d[i][j] = "D"
        # catching edge cases bc sometimes diagonal is not the best move 
        # if the sequences dont match 
        else: 
            if (a[i][j] == a[i][j-1]+ gap):
                d[i][j] = "L"
            elif (a[i][j] == a[i-1][j]+ gap):
                d[i][j] = "U"
            else:
                d[i][j] = "D"
        
# printing to make sure the matrix is okay
#for i in range(len(a)):
#    for j in range(len(a[i])):
#        print((i,j,a[i][j],d[i][j]),end=" ")
#    print(' ')

# these are lists to help with the printing for alignment
seq1_corr = []
seq_mid = []
seq2_corr = []

j = len(seq2)
i = len(seq1)

# this is the starting position
prev_direction = d[i][j]

#score is the bottom right value 
print('score: '+ str(a[i][j]))

# this while loop is for printing the the alignment 
# it looks at the main matrix values and directions and runs through it to print
while(i>0 or j>0):
    # 0 is diagonal 
    #print(i,j)
    if prev_direction == "D": 
        seq1_corr += seq1[i-1]
        if (seq1[i-1] == seq2[j-1]):
            seq_mid += '|'
        else:
            seq_mid += '*'
        seq2_corr += seq2[j-1]
        j = j-1
        i = i-1
    # left
    elif prev_direction == "L" or prev_direction =="A": 
        seq1_corr += '-'
        seq_mid += ' '
        seq2_corr += seq2[j-1]
        j = j - 1
    else: 
        seq1_corr += seq1[i-1]
        seq_mid += ' '
        seq2_corr += '-'
        i = i - 1
    prev_direction = d[i][j] 

# these are the functions for printing the aligned sequences 
seq1_corr.reverse()
seq2_corr.reverse()
seq_mid.reverse()

seq1ToStr = ' '.join([str(elem) for elem in seq1_corr])
seq2ToStr = ' '.join([str(elem) for elem in seq2_corr])
seqMToStr = ' '.join([str(elem) for elem in seq_mid])


print(seq2ToStr)
print(seqMToStr)
print(seq1ToStr)



