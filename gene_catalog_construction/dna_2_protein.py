#!/usr/bin/python

##translate gene catalog

def dna_2_protein(src_filename,faa_filename):
	from Bio import SeqIO
	from Bio.Alphabet import IUPAC
	from Bio.Seq import Seq
	from Bio.Seq import translate
	from Bio.SeqRecord import SeqRecord

	input_handle  = open(src_filename, "r")
	output_handle = open(faa_filename, "w")

	codonTableID = 11

	# Could use this to load all records into memory
	# records = list(SeqIO.parse(src_filename, "fasta"))

	i = 0
	for seq_record in SeqIO.parse(input_handle, "fasta", alphabet=IUPAC.unambiguous_dna):
	    i = i + 1

	    if i% 100 == 0:
	        print("Entry%d:%s"% (i, seq_record.name))

	    # Use this to see help on the translate method:
	    # help (records[0].seq.translate)
	    # You could specify the translation table like this:
	    # seq_record.seq.translate(table="Bacterial")
	    # Note that using cds="true" instructs the code to verify that each sequence is a true CDS

	    try:
	        proteinRecord = SeqRecord(seq_record.seq.translate(cds=False, to_stop=True, table=codonTableID), seq_record.name)
	        proteinRecord.description = seq_record.description

	        SeqIO.write(proteinRecord, output_handle, "fasta")

	        #output_handle.write("&gt;%s%s\n%s\n"% (
	        #   seq_record.name,
	        #   seq_record.description,
	        #   seq_record.seq.translate(cds="false", to_stop="false", table=codonTableID)))

	    except Exception as inst:
	        print("Error translating%s,%s"% (seq_record.name, inst.args[0]))

	        proteinRecord = SeqRecord(translate(seq_record.seq, codonTableID), seq_record.name)
	        proteinRecord.description = seq_record.description + " (translation warning: " + inst.args[0] + ")"

	        SeqIO.write(proteinRecord, output_handle, "fasta")

	        pass

	output_handle.close()
	input_handle.close()
