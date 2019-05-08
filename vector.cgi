#!/usr/local/bin/python3
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline as blastn
from tempfile import NamedTemporaryFile
from Bio.Blast import NCBIWWW
import cgi, cgitb
import jinja2
import re
#import mysql.connector
import os, sys, io

cgitb.enable()



BLAST_EXE = '/var/www/html/zkabir1/finalproject/ncbi-blast-2.6.0+/bin/blastn'  # NETWORK SERVER BLAST FOLDER
#BLAST_EXE = '/home/mkhassan/Downloads/ncbi-blast-2.9.0+/bin/blastn'  # BLAST LOCAL MACHINE

#DB_BASE_PATH = '/home/mkhassan/Downloads/ncbi-blast-2.9.0+/algaedb'  # DATABASE LOCAL LOCATION

DB_BASE_PATH =  '/var/www/html/zkabir1/finalproject/dbases/algaedb'  # custom database server     

MASK = 'N'

def create_rel(XMLin):
    """
    Create a dictionary that relate the sequence name with the region to mask
    """

    bat1 = {}
    output = io.StringIO()
    output.write(XMLin)
    output.seek(0)
    b_records = NCBIXML.parse(output)
    for b_record in b_records:
        for alin in b_record.alignments:
            for hsp in alin.hsps:
                qs = hsp.query_start
                qe = hsp.query_end
                if qs > qe:
                    qe, qs = qs, qe
                record_id = b_record.query.split(' ')[0]
                if record_id not in bat1:
                    bat1[record_id] = [(qs, qe)]
                else:
                    bat1[record_id].append((qs,qe))
    return bat1

def maskseqs(ffh, bat1):
    """
    Take a FASTA file and apply the mask using the positions in the dictionary
    """
    outseqs = []
    for record in SeqIO.parse(ffh, 'fasta'):
        if record.id in bat1:
            #Generate a mutable sequence object to store
            #the sequence with the "mask".
            mutable_seq = record.seq.tomutable()
            coords = bat1[record.id]
            for x in coords:
                mutable_seq[x[0]:x[1]] = MASK*(x[1] - x[0])
            seq_rec = SeqRecord(mutable_seq, record.id,'','')
            outseqs.append(seq_rec)
        else:
            # Leave the sequence as found
            outseqs.append(record)
    return outseqs

#CGI STUFF
# read CGI form values
form = cgi.FieldStorage()
seqs= form.getvalue("seqs")

blastdb = form.getvalue('vector', 'customdb')

#seqs = '>233\nAACGAACGCTGGCGGCAGGCTTAACACATGCAAGTCGAGCGCCCCGCAAGGGGAGCGGCAGACGGGTGAGTAACGCGTGGGAATCTACCGTGCCCTACGGTTGGGCCGTGTCTCAGTCCCAATGTGGCTGATCATCCTCTCAGACCAGCTATGGATCGTCGCCTTGGTAGGCCTTTACCCCACCAACTAGCTAATCCAAC'


#JINJA2 STUFF
# This line tells the template loader where to search for template files
templateLoader = jinja2.FileSystemLoader( searchpath="./templates" )

# This creates your environment and loads a specific template
env = jinja2.Environment(loader=templateLoader)
template = env.get_template('results.html')





#FOR CGI SCRIPT TO USE FILE UPLOAD
# Check if the textarea is empty
#if not seqs:
# Since the textarea is empty, check the uploaded file
#seqs = form.getvalue("seqsfile")

# algae study custom database
#db = os.path.join(DB_BASE_PATH, 'algaedb1')

if blastdb == 'customdb':
	db = os.path.join(DB_BASE_PATH, 'algaedb1')


#if blast_dbase == 'customdb':
    #db = os.path.join(DB_BASE_PATH, 'algaedb1')
#elif blast_dbase == 'ncbivector':
    #db = os.path.join(DB_BASE_PATH, 'algaedb1')
"""
function to connect to servert and do a blast search of sequence on server
def blastweb(seq):
my_query = SeqIO.read("test.fasta", format="fasta")
result_handle = NCBIWWW.qblast("blastn", "nr", my_query.seq)
blast_result = open("my_blast.xml", "w")
blast_result.write(result_handle.read())
blast_result.close()
result_handle.close()
"""
#else:
        #"unexpected string sent, use default customdb"
    #db = os.path.join(DB_BASE_PATH, 'algaedb1')

def result(seqs, db):

    # Create a temporary file
    with NamedTemporaryFile(mode='w') as fasta_in_fh:
        # Write the user entered sequence into this temporary file
        fasta_in_fh.write(seqs)
        # Flush data to disk without closing and deleting the file,
        # since that closing a temporary file also deletes it
        fasta_in_fh.flush()
        # Get the name of the temporary file
        file_in = fasta_in_fh.name
        # Run the BLAST query
        blastn_cline = blastn(cmd=BLAST_EXE, query=file_in, db=db,
                              evalue=.0005, outfmt=5)
        rh, eh = blastn_cline()
        # Create contamination position and store it in a dictionary
        bat1 = create_rel(rh)
        # Get the sequences masked
        newseqs = maskseqs(file_in, bat1)
    with io.StringIO() as fasta_out_fh:
        SeqIO.write(newseqs, fasta_out_fh, 'fasta')
        fasta_out_fh.seek(0)
        finalout = fasta_out_fh.read()
    return finalout
    






print("Content-Type: text/html\n\n")
finalout = result(seqs, db)
print(template.render(finalout=finalout,))
