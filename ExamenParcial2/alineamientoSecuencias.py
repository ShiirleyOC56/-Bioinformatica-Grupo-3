import numpy as np
from Bio import SeqIO
from itertools import product
import os, glob
from math import log, sqrt
os.environ['QT_QPA_PLATFORM']='offscreen'
from ete3 import Tree
from IPython.display import display, Image
from math import log

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

"""def needleman_wunsch(seq1,seq2,Ss=False,match=0,mismatch=0,gap=0):

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
  #return '\n'.join([alineamiento1, alineamiento2])"""
def needleman_wunsch1(seq1,seq2,match,mismatch,gap,S=False):
  len_seq1=len(seq1)
  len_seq2=len(seq2)
  print("-------EN funcion-----------")
  print(match)
  print(mismatch)
  print(gap)
  #creamos la matriz de ceros
  m_inicial=np.zeros((len_seq1+1,len_seq2+1))
  #Llenamos la primera fila y la primera columna de acuerdo al gap
  m_inicial[:,0] = np.linspace(0,len_seq1*gap,len_seq1 + 1)
  m_inicial[0,:] = np.linspace(0,len_seq2*gap,len_seq2 + 1)

  # Scores temporales
  t = np.zeros(3)
  for i in range(len_seq1):
      for j in range(len_seq2):
          t[0]=m_inicial[i,j]+Similitud(seq1[i],seq2[j],S,match,mismatch)          
          t[1] = m_inicial[i,j+1] + gap
          t[2] = m_inicial[i+1,j] + gap
          tmax = np.max(t)
          m_inicial[i+1,j+1] = tmax

  # Trace through an optimal alignment.
  
  alineamiento1=""
  alineamiento2=""
  i = len_seq1
  j = len_seq2
  score=m_inicial[i][j]
  alineamientos=[("","",i,j,score)]
  dic = {}
  print(m_inicial)
  while alineamientos[0][2]>0 and alineamientos[0][3]>0:
    
    for ali in range(len(alineamientos)):

      tupla=alineamientos[ali]

      alineamiento1=tupla[0]
      alineamiento2=tupla[1]
      i = tupla[2]
      j = tupla[3]
      #print(i,"i")
      #print(j,"j")

      score=m_inicial[i][j]

      #print("score",score)
      scoreDiag=m_inicial[i-1][j-1]
      #print("scoreDiag",scoreDiag)
      scoreUp=m_inicial[i-1][j]
      #print("scoscoreUpre",scoreUp)
      scoreLeft=m_inicial[i][j-1]
      #print("scoreLeft",scoreLeft)

      sipi=False
      minimo=(scoreDiag + Similitud(seq1[i-1],seq2[j-1],S,match,mismatch))
      #print(minimo,"minimodiag")
      #print(str(scoreLeft + gap),"minimoleft")
      #print(str(scoreUp+gap),"minimoup")

      if (scoreLeft + gap)> minimo:
        minimo=(scoreLeft + gap)

      if (scoreUp+gap)>minimo:
        minimo=(scoreUp+gap)
      #print(minimo,"minimoFINAL")


      if minimo==(scoreDiag + Similitud(seq1[i-1],seq2[j-1],S,match,mismatch)):
        nuevo=(seq1[i-1]+alineamiento1,seq2[j-1]+alineamiento2,i-1,j-1,tupla[4])
        #print("DIAGONAL")
        if (sipi==False):
          alineamientos[ali]=nuevo
          sipi=True
        else:
          alineamientos.append(nuevo)
      if minimo==(scoreLeft + gap):
        nuevo=(seq1[i-1]+alineamiento1,"-"+alineamiento2,i-1,j,tupla[4])
        #print("No diagonal")
        if (sipi==False):
          alineamientos[ali]=nuevo
          sipi=True
        else:
          alineamientos.append(nuevo)
        
      if minimo==(scoreUp+gap):
        nuevo=("-"+alineamiento1,seq2[j-1]+alineamiento2,i,j-1,tupla[4])
        #print("No diagonal")
        if (sipi==False):
          alineamientos[ali]=nuevo
          sipi=True
        else:
          alineamientos.append(nuevo)
    
    #dic[score]=(alineamiento1,alineamiento2)
    #else:
    #print("sdfsfsdf")
    #print(alineamiento1[::-1])
    #print(alineamiento2[::-1])
  

  for ali in range(len(alineamientos)):

    tupla=alineamientos[ali]

    alineamiento1=tupla[0]
    alineamiento2=tupla[1]
    i = tupla[2]
    j = tupla[3]
    while i>1:
      alineamiento1=seq1[i-1]+alineamiento1
      alineamiento2="-"+alineamiento2
      i-=1
      #dic[score]=(alineamiento1,alineamiento2)
    while j>1:
      alineamiento1="-"+alineamiento1
      alineamiento2=seq2[j-1]+alineamiento2
      j-=1
      #dic[score]=(alineamiento1,alineamiento2)
    
      #print("m_inicial",m_inicial)
    nuevo=(alineamiento1,alineamiento2,i,j,tupla[4])
    alineamientos[ali]=nuevo
  return alineamientos

# PRACTICA 7
def crear_matriz(n,m):
  #creamos la matriz de ceros
  matriz = np.zeros((n+1,m+1))
  return matriz

def finalize(align1, align2,Ss,identicalMatch,mismatch,gap):
    align1 = align1[::-1]    
    align2 = align2[::-1]      
    i,j = 0,0    
    score = 0
    for i in range(0,len(align1)):
        # 
        if align1[i] == align2[i]:                
            score += S(align1[i], align2[i],Ss,identicalMatch,mismatch)
    
        # 
        elif align1[i] != align2[i] and align1[i] != '-' and align2[i] != '-': 
            score += S(align1[i], align2[i],Ss,identicalMatch,mismatch)
    
        #
        elif align1[i] == "-" or align2[i] == "-":          
            score += gap
            
        dic[score]=(align1,align2)
    
    """print ("Score = ", score)
    print (align1)
    print (align2)"""
    return dic

def sw(seq1,seq2,gapcost,identicalMatch,mismatch,Ss = False):
  lens1=len(seq1)
  lens2=len(seq2)
  #creamos la matriz de ceros
  m_inicial=crear_matriz(lens1,lens2)
  punteros = crear_matriz(lens1,lens2)
  m_scores = m_inicial
  T=np.zeros(4) #T[0]:diagonal \, T[1]:arriba |, T[2]: izquierda <-- y T[3]:0
  for i in range(lens1):
    for j in range(lens2):
      T[0]= m_scores[i][j] + Similitud(seq1[i],seq2[j],Ss,identicalMatch,mismatch)
      T[1]= m_scores[i][j+1] + gapcost
      T[2]= m_scores[i+1][j] + gapcost
      tmax = np.max(T)
      m_scores[i+1][j+1] = tmax
      if m_scores[i+1][j+1] == 0:
          punteros[i+1][j+1] = 0 # 0 fin del camino
      if m_scores[i+1][j+1] == T[2]: 
          punteros[i+1][j+1] = 3 # 3 flecha izquierda
      if m_scores[i+1][j+1] == T[1]:
          punteros[i+1][j+1] = 2 # 2 flecha arriba
      if m_scores[i+1][j+1] == T[0]:
          punteros[i+1][j+1] = 1 # 1 FLecha diagonal
           
  print("Matriz de scores") 
  print(m_scores)
  print("Matriz de Punteros")
  print(punteros)
  max_score = 0
  pos_max =np.zeros(2)
  for k in range(lens1+1):
    for l in range(lens2+1):
      if (m_scores[k][l]>= max_score):
        max_score=m_scores[k][l]
        pos_max[0]=k
        pos_max[1]=l
  
  align1=""
  align2=""
  i=int(pos_max[0])
  j=int(pos_max[1])
  dic = {}
  score=0
  # 1 : DIAGONAL- 2: up - 3: left
  while punteros[i][j] != 0:
    if punteros[i][j] == 1:
        align1 += seq1[i-1]
        align2 += seq2[j-1]
        i -= 1
        j -= 1
    elif punteros[i][j] == 2:
        align1 += '-'
        align2 += seq2[j-1]
        j -= 1
    elif punteros[i][j] == 3:
        align1 += seq1[i-1]
        align2 += '-'
        i -= 1 
    ##dic[score]=(align1,align2)

  #finalize(align1,align2,Ss,identicalMatch,mismatch,gapcost)
  align1 = align1[::-1]    
  align2 = align2[::-1]      
  i,j = 0,0    
  for i in range(0,len(align1)):
      # 
      if align1[i] == align2[i]:                
          score += Similitud(align1[i], align2[i],Ss,identicalMatch,mismatch)
  
      # 
      elif align1[i] != align2[i] and align1[i] != '-' and align2[i] != '-': 
          score += Similitud(align1[i], align2[i],Ss,identicalMatch,mismatch)
  
      #
      elif align1[i] == "-" or align2[i] == "-":
          score += gapcost
          
      dic[score]=(align1,align2)
  
  """print ("Score = ", score)
  print (align1)
  print (align2)"""
  print(dic)
  return dic

#BLAST-
fp = []
def bdtoarray(dir):
  folder_path = dir #'bd/'
  fasta_paths = glob.glob(os.path.join(folder_path, '*.fasta'))
  basededatos=[]
  for fasta_path in fasta_paths:
      fp.append(fasta_path)
      for seq_record in SeqIO.parse(fasta_path, "fasta"):
          #print(seq_record.id)
          basededatos.append(str(seq_record.seq.upper()))
  return basededatos

AMINOACID_LIST = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
BLOSUM62 = {
    ('W', 'F'): 1, ('L', 'R'): -2, ('S', 'P'): -1, ('V', 'T'): 0,
    ('Q', 'Q'): 5, ('N', 'A'): -2, ('Z', 'Y'): -2, ('W', 'R'): -3,
    ('Q', 'A'): -1, ('S', 'D'): 0, ('H', 'H'): 8, ('S', 'H'): -1,
    ('H', 'D'): -1, ('L', 'N'): -3, ('W', 'A'): -3, ('Y', 'M'): -1,
    ('G', 'R'): -2, ('Y', 'I'): -1, ('Y', 'E'): -2, ('B', 'Y'): -3,
    ('Y', 'A'): -2, ('V', 'D'): -3, ('B', 'S'): 0, ('Y', 'Y'): 7,
    ('G', 'N'): 0, ('E', 'C'): -4, ('Y', 'Q'): -1, ('Z', 'Z'): 4,
    ('V', 'A'): 0, ('C', 'C'): 9, ('M', 'R'): -1, ('V', 'E'): -2,
    ('T', 'N'): 0, ('P', 'P'): 7, ('V', 'I'): 3, ('V', 'S'): -2,
    ('Z', 'P'): -1, ('V', 'M'): 1, ('T', 'F'): -2, ('V', 'Q'): -2,
    ('K', 'K'): 5, ('P', 'D'): -1, ('I', 'H'): -3, ('I', 'D'): -3,
    ('T', 'R'): -1, ('P', 'L'): -3, ('K', 'G'): -2, ('M', 'N'): -2,
    ('P', 'H'): -2, ('F', 'Q'): -3, ('Z', 'G'): -2, ('X', 'L'): -1,
    ('T', 'M'): -1, ('Z', 'C'): -3, ('X', 'H'): -1, ('D', 'R'): -2,
    ('B', 'W'): -4, ('X', 'D'): -1, ('Z', 'K'): 1, ('F', 'A'): -2,
    ('Z', 'W'): -3, ('F', 'E'): -3, ('D', 'N'): 1, ('B', 'K'): 0,
    ('X', 'X'): -1, ('F', 'I'): 0, ('B', 'G'): -1, ('X', 'T'): 0,
    ('F', 'M'): 0, ('B', 'C'): -3, ('Z', 'I'): -3, ('Z', 'V'): -2,
    ('S', 'S'): 4, ('L', 'Q'): -2, ('W', 'E'): -3, ('Q', 'R'): 1,
    ('N', 'N'): 6, ('W', 'M'): -1, ('Q', 'C'): -3, ('W', 'I'): -3,
    ('S', 'C'): -1, ('L', 'A'): -1, ('S', 'G'): 0, ('L', 'E'): -3,
    ('W', 'Q'): -2, ('H', 'G'): -2, ('S', 'K'): 0, ('Q', 'N'): 0,
    ('N', 'R'): 0, ('H', 'C'): -3, ('Y', 'N'): -2, ('G', 'Q'): -2,
    ('Y', 'F'): 3, ('C', 'A'): 0, ('V', 'L'): 1, ('G', 'E'): -2,
    ('G', 'A'): 0, ('K', 'R'): 2, ('E', 'D'): 2, ('Y', 'R'): -2,
    ('M', 'Q'): 0, ('T', 'I'): -1, ('C', 'D'): -3, ('V', 'F'): -1,
    ('T', 'A'): 0, ('T', 'P'): -1, ('B', 'P'): -2, ('T', 'E'): -1,
    ('V', 'N'): -3, ('P', 'G'): -2, ('M', 'A'): -1, ('K', 'H'): -1,
    ('V', 'R'): -3, ('P', 'C'): -3, ('M', 'E'): -2, ('K', 'L'): -2,
    ('V', 'V'): 4, ('M', 'I'): 1, ('T', 'Q'): -1, ('I', 'G'): -4,
    ('P', 'K'): -1, ('M', 'M'): 5, ('K', 'D'): -1, ('I', 'C'): -1,
    ('Z', 'D'): 1, ('F', 'R'): -3, ('X', 'K'): -1, ('Q', 'D'): 0,
    ('X', 'G'): -1, ('Z', 'L'): -3, ('X', 'C'): -2, ('Z', 'H'): 0,
    ('B', 'L'): -4, ('B', 'H'): 0, ('F', 'F'): 6, ('X', 'W'): -2,
    ('B', 'D'): 4, ('D', 'A'): -2, ('S', 'L'): -2, ('X', 'S'): 0,
    ('F', 'N'): -3, ('S', 'R'): -1, ('W', 'D'): -4, ('V', 'Y'): -1,
    ('W', 'L'): -2, ('H', 'R'): 0, ('W', 'H'): -2, ('H', 'N'): 1,
    ('W', 'T'): -2, ('T', 'T'): 5, ('S', 'F'): -2, ('W', 'P'): -4,
    ('L', 'D'): -4, ('B', 'I'): -3, ('L', 'H'): -3, ('S', 'N'): 1,
    ('B', 'T'): -1, ('L', 'L'): 4, ('Y', 'K'): -2, ('E', 'Q'): 2,
    ('Y', 'G'): -3, ('Z', 'S'): 0, ('Y', 'C'): -2, ('G', 'D'): -1,
    ('B', 'V'): -3, ('E', 'A'): -1, ('Y', 'W'): 2, ('E', 'E'): 5,
    ('Y', 'S'): -2, ('C', 'N'): -3, ('V', 'C'): -1, ('T', 'H'): -2,
    ('P', 'R'): -2, ('V', 'G'): -3, ('T', 'L'): -1, ('V', 'K'): -2,
    ('K', 'Q'): 1, ('R', 'A'): -1, ('I', 'R'): -3, ('T', 'D'): -1,
    ('P', 'F'): -4, ('I', 'N'): -3, ('K', 'I'): -3, ('M', 'D'): -3,
    ('V', 'W'): -3, ('W', 'W'): 11, ('M', 'H'): -2, ('P', 'N'): -2,
    ('K', 'A'): -1, ('M', 'L'): 2, ('K', 'E'): 1, ('Z', 'E'): 4,
    ('X', 'N'): -1, ('Z', 'A'): -1, ('Z', 'M'): -1, ('X', 'F'): -1,
    ('K', 'C'): -3, ('B', 'Q'): 0, ('X', 'B'): -1, ('B', 'M'): -3,
    ('F', 'C'): -2, ('Z', 'Q'): 3, ('X', 'Z'): -1, ('F', 'G'): -3,
    ('B', 'E'): 1, ('X', 'V'): -1, ('F', 'K'): -3, ('B', 'A'): -2,
    ('X', 'R'): -1, ('D', 'D'): 6, ('W', 'G'): -2, ('Z', 'F'): -3,
    ('S', 'Q'): 0, ('W', 'C'): -2, ('W', 'K'): -3, ('H', 'Q'): 0,
    ('L', 'C'): -1, ('W', 'N'): -4, ('S', 'A'): 1, ('L', 'G'): -4,
    ('W', 'S'): -3, ('S', 'E'): 0, ('H', 'E'): 0, ('S', 'I'): -2,
    ('H', 'A'): -2, ('S', 'M'): -1, ('Y', 'L'): -1, ('Y', 'H'): 2,
    ('Y', 'D'): -3, ('E', 'R'): 0, ('X', 'P'): -2, ('G', 'G'): 6,
    ('G', 'C'): -3, ('E', 'N'): 0, ('Y', 'T'): -2, ('Y', 'P'): -3,
    ('T', 'K'): -1, ('A', 'A'): 4, ('P', 'Q'): -1, ('T', 'C'): -1,
    ('V', 'H'): -3, ('T', 'G'): -2, ('I', 'Q'): -3, ('Z', 'T'): -1,
    ('C', 'R'): -3, ('V', 'P'): -2, ('P', 'E'): -1, ('M', 'C'): -1,
    ('K', 'N'): 0, ('I', 'I'): 4, ('P', 'A'): -1, ('M', 'G'): -3,
    ('T', 'S'): 1, ('I', 'E'): -3, ('P', 'M'): -2, ('M', 'K'): -1,
    ('I', 'A'): -1, ('P', 'I'): -3, ('R', 'R'): 5, ('X', 'M'): -1,
    ('L', 'I'): 2, ('X', 'I'): -1, ('Z', 'B'): 1, ('X', 'E'): -1,
    ('Z', 'N'): 0, ('X', 'A'): 0, ('B', 'R'): -1, ('B', 'N'): 3,
    ('F', 'D'): -3, ('X', 'Y'): -1, ('Z', 'R'): 0, ('F', 'H'): -1,
    ('B', 'F'): -3, ('F', 'L'): 0, ('X', 'Q'): -1, ('B', 'B'): 4
}

def comparar(a,b):
   if (a, b) in BLOSUM62.keys():
        return BLOSUM62[(a, b)]
   else:
        return BLOSUM62[(b, a)]

def score(a,b):
  return comparar(a[0],b[0])+comparar(a[1],b[1])+comparar(a[2],b[2])

def BLAST(cad1,DB):
  dic = {}
  #DB=["PQLPITNFSRDWQSGRALGALVDSCAEYYPMVPDSWDASKPVTNAREAMQQADDWLGIPQ","VITPEEIVDPNVDEHSVMTYLSQFPKAKLKPGAPLRPKLNPKKARAYGPGIEPTGNMVKK","RAEFTVETRSAGQGEVLVYVEDPAGHQEEAKVTANNDKNRTFSVWYVPEVTGTHKVTVLF","AGQHIAKSPFEVYVDKSQGDASKVTAQGPGLEPSGNIANKTTYFEIFTAGAGTGEVEVVI","QDPMGQKGTVEPQLEARGDSTYRCSYQPTMEGVHTVHVTFAGVPIPRSPYTVTVGQACNP","SACRAVGRGLQPKGVRVKETADFKVYTKGAGSGELKVTVKGPKGEERVKQKDLGDGVYGF"]
  lencad1=len(cad1)
  partes=[[cad1[i:i+3]] for i in range(0,lencad1-2)]

  # vECINOS

  posibilidades=[''.join(j) for j in  product('CSTPAGNDEQHRKMILVFYW', repeat=3)]
  #print("PARTES",partes)
  for i in partes:
    for j in posibilidades:
      if ((comparar(i[0][0], j[0])+comparar(i[0][1], j[1])+comparar(i[0][2], j[2]))>12):
        scor=comparar(i[0][0], j[0])+comparar(i[0][1], j[1])+comparar(i[0][2], j[2])
        #print('score=',scor)
        i.append(j)
    i.remove(i[0])
  #print("PARTES",partes)

  #BLAST

  for k in range(len(partes)-1): # Recorre el vector con los vecinos [[abc,aeb,acc..],[],[]]
    for kk in partes[k]: # Recorre el subvector [abc,aeb,acc..]
      cont=0
      for i in DB: # Recorre ["","",""]
        for j in range(len(i)-2): #Recorre "ABCDEFGHI" DESDE "A" HASTA "G" 
          if kk[0]==i[j]:
            if kk[1]==i[j+1]:
              if kk[2]==i[j+2]:
                x = []
                print("----------------")
                print("KK",kk)
                print("CAD1",cad1[k:k+3])
                print("Secuencia ",fp[cont])
                print(i[j:j+3])
                ttt=score(cad1[k:k+3],i[j:j+3])
                print("score",ttt)
                izqQuery=k
                derQuery=k+2
                izqDB=j
                derDB=j+2
                print("iq",izqQuery)
                print("dq",derQuery)
                print("id",izqDB)
                print("dd",derDB)
                while (ttt>10):#22 PARA PROTINAS Y 20 PARA DNA
                  beneficioDer=0
                  beneficioIzq=0

                  if izqQuery>0 and izqDB>0:
                    beneficioIzq=comparar(cad1[izqQuery-1],i[izqDB-1])
                  else:
                    beneficioIzq=-1000
                  if derQuery<len(cad1)-1 and derDB<len(i)-1:
                    #print("d",derecha)
                    #print("kk",cad1)
                    #print("i",i)
                    beneficioDer=comparar(cad1[derQuery+1],i[derDB+1])
                  else:
                    beneficioDer=-1000
                  
                  if (beneficioDer>beneficioIzq):
                    derQuery+=1
                    derDB+=1
                    ttt+=beneficioDer
                  elif (beneficioDer<beneficioIzq):
                    izqQuery-=1
                    izqDB-=1
                    ttt+=beneficioIzq
                  else:
                    break

                print("cad1",cad1[izqQuery:derQuery+1])
                print("cadi",i[izqDB:derDB+1])
                print("scoreFinal",ttt)
                #x=["triplete inicial","secuencia","cad1","cadi"]
                x.append(kk)
                x.append(fp[cont])
                x.append(cad1[izqQuery:derQuery+1])
                x.append(i[izqDB:derDB+1])
                dic[ttt] = x
        cont+=1
  return dic

#MUSCLE 
def MUSCLE(sequences):
    records = [SeqRecord(Seq(sequence)) for sequence in sequences]
    muscle_cline = MuscleCommandline(clwstrict=True, cmd='muscle')

    handle = StringIO()
    SeqIO.write(records, handle, "fasta")
    data = handle.getvalue()

    stdout, stderr = muscle_cline(stdin=data)
    align = AlignIO.read(StringIO(stdout), "clustal")

    length = align.get_alignment_length()

    return [str(entry.seq) for entry in align]

#jukes cantor

import re, math

def distance (seq1, seq2):
    p = percent_difference_of_nucleotides(seq1, seq2)
    return -0.75 * math.log(1 - (p*4/3)) if p else 0


def percent_difference_of_nucleotides (seq1, seq2, nucleobases=set('ACGT')):
	# percentage of nucleotide difference in two sequences

	diff_count = 0 # number of nucleotide differences
	valid_nucleotides_count = 0.0 # number of valid nucleotides (value is float for computing percentage)

	for a, b in zip(seq1, seq2):
		if a in nucleobases and b in nucleobases:
			valid_nucleotides_count += 1
			if a != b: diff_count += 1
	
	return diff_count / valid_nucleotides_count if valid_nucleotides_count else 0


def jukes_cantor(dir):
  resultados =[]
  pattern = re.compile(r'(?<=\>)(.+?)(?=\|)') # regex to find sequence id
  sequences_read = {} # store previously read sequences  resultados = []
  with open(dir, 'r') as input_file:
    # read fasta data
    while True: 
      line = input_file.readline().rstrip()
      if not line: break
      
      # extract sequence id and read the sequence
      new_seq_id = pattern.search(line).group()
      new_seq = input_file.readline().rstrip()
      
      # compute distance with regard to the previously read sequences
      
      for seq_id in sequences_read:
        item =[]
        d = distance(new_seq, sequences_read[seq_id])

        # discard distances greater than 4%
        if d <= .04:
            item.append(new_seq_id)
            item.append(seq_id)
            item.append(d)
            print(f'{new_seq_id}\t{seq_id}\t{d}')
            resultados.append(item)
     
      # add new sequence to previously read sequences
      sequences_read.update({new_seq_id: new_seq})
    return resultados

#KIMURA
def K2Pdistance(seq1,seq2):
    """
    Kimura 2-Parameter distance = -0.5 log( (1 - 2p -q) * sqrt( 1 - 2q ) )
    where:
    p = transition frequency
    q = transversion frequency
    """
    pairs = []

    #collect ungapped pairs
    for x in zip(seq1,seq2):
        if '-' not in x: pairs.append(x)
        
    ts_count=0
    tv_count=0
    length = len(pairs)
    
    transitions = [ "AG", "GA", "CT", "TC"]
    transversions = [ "AC", "CA", "AT", "TA",
                      "GC", "CG", "GT", "TG" ]

    for (x,y) in pairs:
        if x+y in transitions: ts_count += 1 
        elif x+y in transversions: tv_count += 1
    
    p = float(ts_count) / length
    q = float(tv_count) / length
    try: d = -0.5 * log( (1 - 2*p - q) * sqrt( 1 - 2*q ) )
    except ValueError: 
        return None
    return d 

#TAMURA DISTANCE
def Tamuradistance(seq1,seq2):
    pairs = []
    
    #collect ungapped pairs
    for x in zip(seq1,seq2):
        if '-' not in x: pairs.append(x)
        
    ts_count=0
    tv_count=0
    length = len(pairs)
    
    transitions = [ "AG", "GA", "CT", "TC"]
    transversions = [ "AC", "CA", "AT", "TA",
                      "GC", "CG", "GT", "TG" ]

    for (x,y) in pairs:
        if x+y in transitions: ts_count += 1 
        elif x+y in transversions: tv_count += 1
    
    p = float(ts_count) / length
    q = float(tv_count) / length
    gc1 = sum(estimate_nucleotide_frequencies(seq1)[1:3])
    gc2 = sum(estimate_nucleotide_frequencies(seq2)[1:3])
    c = gc1 + gc2 - 2 * gc1 * gc2

    try: d = -c * log( 1 - p/c - q) - 0.5 * ( 1 - c ) * log ( 1 - 2*q )
    except ValueError: 
        return None
    return d

#TAJIMA 
def estimate_nucleotide_frequencies(seq):
    seq = seq.replace('-','').upper()
    A = seq.count('A')
    C = seq.count('C')
    G = seq.count('G')
    T = seq.count('T')
    length = float(len(seq))
    return [ x/length for x in [A,C,G,T] ]

def pdistance(seq1, seq2):
    p = 0
    pairs = []
    for x in zip(seq1,seq2):
        if '-' not in x: pairs.append(x)
    #for (x,y) in zip(seq1,seq2):
    for (x,y) in pairs:
        if x != y:
            p += 1
    #length = (len(seq1) + len(seq2)) / 2
    length = len(pairs)
    return float(p) / length

def TNdistance(seq1, seq2):
    ns = ['A','C','G','T']
    G = estimate_nucleotide_frequencies(seq1 + seq2)
    p = pdistance(seq1,seq2)
    pairs = []
    h = 0

    #collect ungapped pairs
    for x in zip(seq1,seq2):
        if '-' not in x: pairs.append(x)
       
    #pair frequencies are calculated for AC, AG, AT, CG, CT, GT (and reverse order)
    for i in range(len(ns)-1):
        for j in range(i+1,len(ns)):
            if i != j: paircount = pairs.count( (ns[i], ns[j]) ) + pairs.count( (ns[j], ns[i]) )
            Xij_sq = (float(paircount)/len(pairs))**2
            GiGj = G[i]*G[j]
            h += 0.5*Xij_sq/GiGj  #h value used to calculate b
    
    b = 0.5*(1-sum([x**2 for x in G])+p**2/h)
    try: d = -b * log(1 - p/b)
    except ValueError: 
        return None
    return d

#UPGMA
def lowest_cell(table):
    
    min_cell = float("inf")
    x, y = -1, -1

    for i in range(len(table)):
        for j in range(len(table[i])):
            if table[i][j] < min_cell:
                min_cell = table[i][j]
                x, y ,w= i, j,min_cell
    return x, y,w



def join_labels(labels, a, b,w):
    
    if b < a:
        a, b = b, a
    
    labels[a]=labels[a]+":"+str(w/2)
    labels[b]=labels[b]+":"+str(w/2)

    labels[a] = "(" + labels[a] + "," + labels[b]+ " )"

    del labels[b]


def join_table(table, a, b):

    if b < a:
        a, b = b, a

    row = []
    for i in range(0, a):
        row.append((table[a][i] + table[b][i])/2)
    table[a] = row

    for i in range(a+1, b):
        table[i][a] = (table[i][a]+table[b][i])/2

    for i in range(b+1, len(table)):
        table[i][a] = (table[i][a]+table[i][b])/2

        del table[i][b]
    del table[b]



def UPGMA(table, labels):

    while len(labels) > 1:

        x, y,w = lowest_cell(table)

        join_table(table, x, y)
        join_labels(labels, x, y,w)
        #print(table)

    return labels[0]



def alpha_labels(start, end):
    labels = []
    for i in range(ord(start), ord(end)+1):
        labels.append(chr(i))
    return labels

M_labels = alpha_labels("A", "D")   
M = [[],[8],
        [4, 8],
        [6, 8, 6]]

print(UPGMA(M, M_labels)) 


#NEGIBOR JOINING
def neighborjoining(dis_map, n, L):
	if n==2:
		print (dis_map)
		return dis_map
	else:
		#calculate the r_i coefficients
		r = {}
		for (i,j) in dis_map:
			r[i] = 0

		for (i,j) in dis_map:
			r[i] += float(1/(float(n-2)))*float(dis_map[(i,j)])

		#calculate the Dij coefficients
		D = {}
		for (i,j) in dis_map:
				D[(i,j)] = dis_map[(i,j)]-r[i]-r[j]

		min_D = float('inf')
		min_i = -1
		min_j = -1
		for (i,j) in D:
			val = D[(i,j)]

			if val < min_D and i!=j:
				min_D = val
				min_i = i
				min_j = j

		print (min_i+1, min_j+1, L+1)

		new_dis_map = {}
		k = L
		for (i,j) in dis_map:
			#print min_i,min_j, i ,j
			if j==min_i or j==min_j:
				pass
			else:
				new_dis_map[(k,j)] = .5 * float(dis_map[(min_i,j)]+dis_map[(min_j,j)]-dis_map[(min_i,min_j)])
				new_dis_map[(j,k)] = new_dis_map[(k,j)]
		new_dis_map[(k,k)] = 0

		#remove min_i and min_j from dismap
		for (i,j) in dis_map:
			if i == min_i or j==min_j or i==min_j or j==min_i:
				pass
			else:
				#print i,j,min_i,min_j
				new_dis_map[(i,j)] = dis_map[(i,j)]

		#print new_dis_map
		d_ik = .5*float(dis_map[(min_i,min_j)] + r[min_i] - r[min_j])
		d_jk = dis_map[(min_i,min_j)] - d_ik

		print (d_ik, d_jk)

		#print new_dis_map
		neighborjoining(new_dis_map,n-1,L+1)

"""dis_data = open('data.txt')

dis_map = {}
i = 0
for line in dis_data:
	split_line = line.split()
	print (split_line)
	for j in range(len(split_line)):
		dis_map[(i,j)] = int(split_line[j])
	i+=1
dis_data.close()

neighborjoining(dis_map, 8,8)"""