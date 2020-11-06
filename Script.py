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

def concatenate_and_get_reverse_of_complement(secuencia1, secuencia2):
    secuencia1 = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
    secuencia2 = Seq("GACGTACCTATGTATAGCGATACGTTAGCTAC")
    secuencia = secuencia1 + secuencia2
    print(secuencia.reverse_complement())
concatenate_and_get_reverse_of_complement(secuencia1, secuencia2)
