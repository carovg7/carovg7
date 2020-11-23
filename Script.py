from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

# Se debe de cambiar a la dirección del nuevo archivo .gbk o .fasta a leer
filename = "/home/carolina/carovg7/data/m_cold.fasta"

def summarize_contents(filename):
	lista = []
	lista = os.path.split(filename)
	cadena = " "
	cadena = ("\nfile: "+ lista[1] + "\npath: " + lista[0])
	File_Extension = []
	File_Extension = os.path.splitext(filename)
	if(File_Extension[1]==".gbk"):
		type_file="genbank"
	else:
		type_file="fasta"
		pass
	records = list(SeqIO.parse(filename, type_file))
	cadena += ("\nnum_records: " + str(len(records)))
	cadena += ("\nrecord(s): ")
	for seq_record in SeqIO.parse(filename, type_file):
		cadena += ("\n- id: " + str(seq_record.id))
		cadena += ("\nname: " + seq_record.name)
		cadena += ("\ndescription: " + str(seq_record.description))
	return cadena
if __name__=="__main__":
	resultado = summarize_contents(filename)
	print(resultado)

#///////////////////////////////////////////////////////////////////////
from Bio.Seq import Seq

secuencia1 = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
secuencia2 = Seq("GACGTACCTATGTATAGCGATACGTTAGCTAC")

def concatenate_and_get_reverse_of_complement(secuencia1, secuencia2):
    secuencia = secuencia1 + secuencia2
    rev_com = secuencia.reverse_complement()
    return rev_com
    print(rev_com)
concatenate_and_get_reverse_of_complement(secuencia1, secuencia2)

#/////////////////////////////////////////////////////////////////////
from Bio.Seq import Seq

cadena = "ATGTCACTTACTTACTCACAGTCTTAA"
def print_protein_and_stop_codon_using_standard_table(cadena):
    secuencia = Seq(cadena)
    diccionario = {}
    mrna = secuencia.transcribe()
    diccionario['mRNA'] = mrna
    for i in range(len(secuencia)):
            if((secuencia[i*3:i*3+3] == 'ATG') or (secuencia[i*3:i*3+3] == 'TTG') or (secuencia[i*3:i*3+3] == 'CTG')):
                proteins = secuencia[i*3:].translate(to_stop = True)
                diccionario['Proteins'] = proteins

                for j in range(len(secuencia)):
                    if((secuencia[j*3:j*3+3] == 'TAG') or (secuencia[j*3:j*3+3] == 'TAA') or (secuencia[j*3:j*3+3] == 'TGA')):
                        diccionario['Stop codon'] = secuencia[j*3:j*3+3]
                        break

            if(i+1 == len(secuencia)):
                break
    return diccionario
r = print_protein_and_stop_codon_using_standard_table(cadena)
print(r)

#////////////////////////////////////////////////////////////////
from Bio.Data import CodonTable

cadena = "ATATCCACTTAA"
def print_proteins_and_codons_using_mitocondrial_yeast_table(cadena):
    secuencia = Seq(cadena)
    diccionario = {}
    mrna = secuencia.transcribe()
    diccionario['mRNA'] = mrna
    for i in range(len(secuencia)):
            if((secuencia[i*3:i*3+3] == 'ATA') or (secuencia[i*3:i*3+3] == 'ATG') or (secuencia[i*3:i*3+3] == 'GTG')):
                proteins = secuencia[i*3:].translate(table = 3, to_stop = True)
                diccionario['Proteins'] = proteins

                for j in range(len(secuencia)):
                    if((secuencia[j*3:j*3+3] == 'TAA') or (secuencia[j*3:j*3+3] == 'TAG')):
                        diccionario['Stop codon'] = secuencia[j*3:j*3+3]
                        break

            if(i+1 == len(secuencia)):
                break
    return diccionario
r = print_proteins_and_codons_using_mitocondrial_yeast_table(cadena)
print(r)
#///////////////////////////////////////////////////////////////////////////
#FUNCIÓN QUE DADO UN ARCHIVO FASTA CON MÚLTIPLES SECUENCIAS, IMPRIME UN ARCHIVO FASTA PARA CADA RECORD
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

filename = "/home/carolina/data/sequences.fasta"
def extract_sequences(filename):
    cont_files = 0
    param = list(SeqIO.parse(filename, "fasta"))
    for i, record in enumerate(param):
        c = open(f'sequences{i}.fasta', 'w')
        c.write(f'>{record.id}\n')
        c.write(f'{record.seq}')
        c.close()
        cont_files = cont_files + 1
    return cont_files
extract_sequences(filename)
#//////////////////////////////////////////////////////////////////////////
#FUNCIÓN PARA EL REVERSO COMPLEMENTARIO, DEVUELVE UN ARCHIVO

filename = "/home/carolina/data/sequences.fasta"
def extract_sequences_revcomp(filename):
    param = list(SeqIO.parse(filename, "fasta"))
    c = open(f'sequences_revcomp.fasta', 'w')
    for i, record in enumerate(param):
        c.write(f'>{record.id}\n')
        c.write(f'{record.seq.reverse_complement()}\n')
    c.close()
    numrecords = len(param)
    return numrecords
extract_sequences_revcomp(filename)
#/////////////////////////////////////////////////////////////////////////
#FUNCIÓN QUE DADO UN FILE FASTA, IMPRIME SUS RECORDS POR SEPARADO EN FORMATO GBK
filename = "/home/carolina/data/sequences.fasta"
salida = "genbank"
def extract_sequences(filename, salida):
    SeqIO.convert(filename, "fasta", "sequences.gbk", "genbank", molecule_type = "protein")
    param = list(SeqIO.parse("sequences.gbk", "genbank"))
    cont_files = 0
    for i, record in enumerate(param):
        c = open(f"sequences{i}.gbk", "w")
        c.write(record.format("genbank"))
        c.close()
        cont_files = cont_files + 1
    return cont_files
extract_sequences(filename, salida)
#////////////////////////////////////////////////////////////////////////////
