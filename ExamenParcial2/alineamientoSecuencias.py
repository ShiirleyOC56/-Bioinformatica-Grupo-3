import numpy as np
from Bio import SeqIO

def fastatoString(archivoFasta):
  sequences = SeqIO.parse(archivoFasta, "fasta")
  for record in sequences:
    data1 = str(record.seq.upper()) # the fasta file just have one sequence
  return data1  

##ALINEAMIENTO GLOBAL - P6
def Similitud(a,b,S,identicalMatch,mismatch):
  if (S == True):
    if (a == b):
      return 2
    elif ((a=="A" and b=="G") or (a=="C" and b=="T") or (a=="G" and b=="A") or (a=="T" and b=="C")):
      return -5
    else:
      return -7
  else:
    match = identicalMatch
    difmatch = mismatch
    if (a == b):
      return match
    else:
      return difmatch

def needleman_wunsch(seq1,seq2,Ss=False,match=0,mismatch=0,gap=0):

  len_seq1=len(seq1)
  len_seq2=len(seq2)
  #creamos la matriz de ceros
  m_inicial=np.zeros((len_seq1+1,len_seq2+1))
  #Llenamos la primera fila y la primera columna de acuerdo al gap
  m_inicial[:,0] = np.linspace(0,len_seq1*gap,len_seq1 + 1)
  m_inicial[0,:] = np.linspace(0,len_seq2*gap,len_seq2 + 1)

  # Scores temporales
  t = np.zeros(3)
  for i in range(len_seq1):
      for j in range(len_seq2):
          t[0]=m_inicial[i,j]+Similitud(seq1[i],seq2[j],Ss,match,mismatch)          
          t[1] = m_inicial[i,j+1] + gap
          t[2] = m_inicial[i+1,j] + gap
          tmax = np.max(t)
          m_inicial[i+1,j+1] = tmax

  # Trace through an optimal alignment.
  
  alineamiento1=""
  alineamiento2=""
  i = len_seq1
  j = len_seq2
  score=0
  dic = {}

  while i>0 and j>0:
    score=m_inicial[i][j]

    print("score",score)
    scoreDiag=m_inicial[i-1][j-1]
    print("scoreDiag",scoreDiag)
    scoreUp=m_inicial[i][j-1]
    print("scoscoreUpre",scoreUp)
    scoreLeft=m_inicial[i-1][j]
    print("scoreLeft",scoreLeft)

    if score==(scoreDiag + Similitud(seq1[i-1],seq2[j-1],Ss,match,mismatch)):
      alineamiento1=seq1[i-1]+alineamiento1
      alineamiento2=seq2[j-1]+alineamiento2
      i-=1
      j-=1
    elif score==(scoreLeft + gap):
      alineamiento1=seq1[i-1]+alineamiento1
      alineamiento2="-"+alineamiento2
      i-=1
    elif score==(scoreUp+gap):
      alineamiento1="-"+alineamiento1
      alineamiento2=seq2[j-1]+alineamiento2
      j-=1
    
    dic[score]=(alineamiento1,alineamiento2)
    #else:
    #print("sdfsfsdf")
    print(alineamiento1[::-1])
    print(alineamiento2[::-1])

  while i>0:
    alineamiento1=seq1[i-1]+alineamiento1
    alineamiento2="-"+alineamiento2
    i-=1
    dic[score]=(alineamiento1,alineamiento2)
  while j>0:
    alineamiento1="-"+alineamiento1
    alineamiento2=seq2[j-1]+alineamiento2
    j-=1
    dic[score]=(alineamiento1,alineamiento2)
  
  return dic
  #print(dic)
  #print("m_inicial",m_inicial)
  #return '\n'.join([alineamiento1, alineamiento2])