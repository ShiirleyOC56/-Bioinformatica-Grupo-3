# dot_matrix.py

from Bio import SeqIO
import matplotlib.pyplot as plt
sequences = SeqIO.parse("M3YD12.fasta", "fasta")
for record in sequences:
    data1 = str(record.seq.upper()) # the fasta file just have one sequence
sequences = SeqIO.parse("G3QNL9.fasta", "fasta")
for record in sequences:
    data2 = str(record.seq.upper()) # the fasta file just have one sequence


def sacar_modulo(num, modulo):
  n = num % modulo
  if (n==0):
    return n
  else:
    return n

def completar_cadenas(cadena,ws):
  lista = list(cadena)
  tam=len(lista)  
  res = ws - (sacar_modulo(tam,ws))
  for i in range(res):
    lista.append("")
  
  return lista

def agrupar_window_size(ws,data):
  dataN=completar_cadenas(data,ws)
  dataR=[]
  for i in range(0,len(dataN),ws):
    sublista=[]
    for j in range(ws):
      sublista.append(dataN[i+j])
    dataR.append(sublista)
  return dataR
    
def comparar(cad1,cad2):
  similares=0
  for i in range(len(cad1)):
    if(cad1[i]==cad2[i]):
      similares=similares+1

  return similares

def sacar_porcentaje(similares,total):
  porcentaje=(similares*100)/total
  return porcentaje

def dotplot(data1,data2,ws):
  seq1 = agrupar_window_size(ws,data1)
  seq2 = agrupar_window_size(ws,data2)
  for i in range(len(seq1)):
    for j in range(len(seq2)):
        sim = comparar(seq1[i],seq2[j])
        porce=sacar_porcentaje(sim,ws)
        #TRESHOLD 100 %
        if(porce==100):
          plt.plot(i,j,'g.',label = "Treeshol 100%" , linewidth=2)
        #TRESHOLD 80 %
        if(porce==80):
          plt.plot(i,j,'b|',label = "Treeshol 80%" , linewidth=4)
        #TRESHOLD 60 %
        if(porce==60):
          plt.plot(i,j,'r*',label = "Treeshol 60%" , linewidth=2)
        #TRESHOLD 20 %
        if(porce==20):
          plt.plot(i,j,'y.',label = "Treeshol 20%" , linewidth=2)

  plt.show()
    

dotplot(data1,data2,10)

