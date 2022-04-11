# from Bio import SeqIO
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

# if using biopython this is the code
#record1 = list(SeqIO.parse("seq2.fa", "fasta"))
#seq1 = list(record1[0].seq)
#print(seq1)

#r2 = list(SeqIO.parse("seq1.fa", "fasta"))
#seq2 = list(r2[0].seq)


negative_infinity = float('-inf')
#created matrix a and initalized all values to 0 
M = [[negative_infinity for x in range(len(seq2)+1)] for x in range(len(seq1)+1)]
X = [[negative_infinity for x in range(len(seq2)+1)] for x in range(len(seq1)+1)]
Y = [[negative_infinity for x in range(len(seq2)+1)] for x in range(len(seq1)+1)]
DM = [[("graph") for x in range(len(seq2)+1)] for x in range(len(seq1)+1)]
DX = [[("graph") for x in range(len(seq2)+1)] for x in range(len(seq1)+1)]
DY = [[("graph") for x in range(len(seq2)+1)] for x in range(len(seq1)+1)]

#gap penalites 
h = -5
g = -1

#initalizing 
M[0][0] = 0
X[0][0] = h
Y[0][0] = h

def match (i,j):
    #print((seq1[i-1],seq2[j-1]))
    if (seq1[i-1] == seq2[j-1]): 
        return(1)
    else: 
        return(-2)

# this matrix works in accordance with the algorithm presented in class 
# DM is the direction matrix which states where the maximum value was found 
for i in range(len(M)):
    for j in range(len(M[i])):
        # do somethings for M to understand what direction to move 
        if(i>0 and j>0):
            M[i][j] = max(M[i-1][j-1] + match(i,j), X[i-1][j-1] + match(i,j), Y[i-1][j-1] + match(i,j))

            if (M[i][j] == M[i-1][j-1] + match(i,j)):
                DM[i][j] = "M"
            elif (M[i][j] == X[i-1][j-1] + match(i,j)):
                DM[i][j] = "X"
            else:
                DM[i][j] = "Y"

        # do something for X
        if(i>0):
            X[i][j] = max(M[i-1][j] + h + g, X[i-1][j] + g)

            if (X[i][j] == M[i-1][j] + h+g):
                DX[i][j] = "M"
            else:
                DX[i][j] = "X"


        if(j>0):
            # do something for Y
            Y[i][j] = max(M[i][j-1] + h + g, Y[i][j-1] + g)

            if (Y[i][j] == M[i][j-1] + h+g):
                DY[i][j] = "M"
            else:
                DY[i][j] = "Y"
        

''' print("M:")
for i in range(len(M)):
    for j in range(len(M[i])):
        print((M[i][j],DM[i][j]), end=" ")
    print()

print("X:")
for i in range(len(M)):
    for j in range(len(M[i])):
        print((X[i][j],DX[i][j]), end=" ")
    print()

print("Y:")
for i in range(len(M)):
    for j in range(len(M[i])):
        print((Y[i][j],DY[i][j]), end=" ")
    print() '''

seq1_corr = []
seq_mid = []
seq2_corr = []

#starting positions
j = len(seq2)
i = len(seq1)

#find the highest value of the three arrays to find where to start 
start = max (M[i][j], X[i][j], Y[i][j])

if(start == M[i][j]): 
    direction = DM[i][j]
elif(start == X[i][j]): 
    direction = DX[i][j]
elif(start == Y[i][j]): 
    direction = DY[i][j]

# use the previous direction to understand where the value came from 
# this is useful for printing values 
prev_direction = direction
print('score: ' +str(start))
while(i>0 and j>0):
    if(prev_direction =="X"): 
        if (i==1):
            seq2_corr += seq2[0]
            seq_mid += '|'
        # print(prev_direction, M[i][j])
        else: 
            seq_mid += ' '
            seq2_corr += '-'
        seq1_corr += seq1[i-1]
        i = i - 1
    # if direction is Y, then gap in the seq2 
    elif(prev_direction =="Y"): 
        # print(prev_direction, M[i][j])
        if (j==1):
            seq1_corr += seq1[0]
            seq_mid += '|'
        else:
            seq1_corr += '-'
            seq_mid += ' '
        seq2_corr += seq2[j-1]
        j = j - 1
    #diagonal move 
    else: 
        #print(prev_direction, M[i][j])
        seq1_corr += seq1[i-1]
        seq2_corr += seq2[j-1]
        if (seq1[i-1] == seq2[j-1]):
            seq_mid += '|'
        else:
            seq_mid += '*'
        j = j-1
        i = i-1
    prev_direction = direction
    direction = DM[i][j]
    
seq1_corr.reverse()
seq2_corr.reverse()
seq_mid.reverse()

seq1ToStr = ' '.join([str(elem) for elem in seq1_corr])
seq2ToStr = ' '.join([str(elem) for elem in seq2_corr])
seqMToStr = ' '.join([str(elem) for elem in seq_mid])

print(seq2ToStr)
print(seqMToStr)
print(seq1ToStr)