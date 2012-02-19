import re
import numpy
import subprocess
from os import path, makedirs, listdir
from Bio import SeqIO, GenBank
from Bio.Blast import NCBIWWW
from Bio.Alphabet import generic_dna
from Bio.Blast.Applications import NcbiblastnCommandline, \
    NcbirpsblastCommandline, NcbiblastpCommandline, NcbitblastxCommandline, \
    NcbiblastxCommandline, NcbitblastnCommandline

def ensure_dir(dir_list):
    """Check that the directory exists; if not, create it."""
    for dir_path in dir_list:
        abs_path = path.abspath(dir_path)
        if not path.exists(abs_path):
            try: makedirs(abs_path)
            except Exception as message:
                status = 1
                # TODO: make graceful fail or request input if interactive mode
            else:
                message = 'created path'
                status = 0
        else:
            message = 'path exists'
            status = 0
        report = {'message': message, 'status': status}

def load_multifasta(seqfile):
    """Load multiple records from Fasta file."""
    input_handle = open(seqfile, 'rU')
    multifasta = SeqIO.parse(input_handle, 'fasta')
    fasta_list = list(multifasta)
    for record in fasta_list:
        assert record.id
    return fasta_list

def load_fasta(seqfile):
    """Load single-record Fasta file."""
    input_handle = open(seqfile, 'rU')
    fasta_record = SeqIO.read(input_handle, 'fasta')
    assert fasta_record.id
    input_handle.close()
    return fasta_record

def load_genbank(seqfile):
    """Load single-record GenBank file."""
    parser = GenBank.FeatureParser()
    input_handle = open(seqfile, 'rU')
    gb_record = parser.parse(input_handle)
    input_handle.close()
    return gb_record

def from_dir(ori_dir, pattern):
    """Load filenames in a directory using a regex."""
    # get directory contents
    contents = listdir(ori_dir)
    filenames = []
    for item in contents:
        match = re.match(pattern, item)
        if match:
            filenames.append(item)
    return filenames

def write_genbank(filename, seqrecords):
    """Write GenBank file."""
    output_handle = open(filename, 'w')
    counter = SeqIO.write(seqrecords, output_handle, 'genbank')
    output_handle.close()
    return counter
    
def write_fasta(filename, seqrecords):
    """Write Fasta file."""
    output_handle = open(filename, 'w')
    count = SeqIO.write(seqrecords, output_handle, 'fasta')
    output_handle.close()
    return count

def fas2gbk(fas_file):
    """Convert a FastA file to Genbank format."""
    record = load_fasta(fas_file)
    gbk_file = fas_file[:fas_file.find('.fas')]+'.gbk'
#    record.name = rec_name
#    record.id = rec_name
    record.seq.alphabet = generic_dna
    write_genbank(gbk_file, record)
    return gbk_file

def gbk2fas(gbk_file):
    """Convert a FastA file to Genbank format."""
    record = load_genbank(gbk_file)
    fas_file = gbk_file[:gbk_file.find('.gbk')]+'.fas'
#    record.name = rec_name
#    record.id = rec_name
    write_fasta(gbk_file, record)
    return fas_file

def train_prodigal(seq_file, training_file, mode):
    """Train Prodigal on the entire dataset."""
    cline = "prodigal "+mode+" -i "+seq_file+" -t "+training_file
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    return output

def run_prodigal(in_file, an_gbk, an_aa, trn_file, mode):
    """Annotate sequence records individually using Prodigal."""
    cline = "prodigal "+mode+" -i "+in_file\
                            +" -o "+an_gbk\
                            +" -a "+an_aa\
                            +" -t "+trn_file
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    return output

def make_blastDB(name, infile, db_type):
    """Make BLAST database from FASTA input file."""
    cline = "makeblastdb -in "+ infile +" -dbtype "+ db_type +" -title "+  \
            infile +" -out "+ name +" -parse_seqids"
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    return output

def local_blastn_2file(query_file, dbfile_path, outfile, prefs):
    """Perform blastn against local database."""
    cline = NcbiblastnCommandline(query=query_file,
                                  db=dbfile_path,
                                  out=outfile,
                                  evalue=prefs['evalue'],
                                  outfmt=prefs['outfmt_pref'])
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate() # forces the main script to wait

def local_blastp_2file(query_file, dbfile_path, outfile, prefs):
    """Perform blastp against local database."""
    cline = NcbiblastpCommandline(query=query_file,
                                  db=dbfile_path,
                                  out=outfile,
                                  evalue=prefs['evalue'],
                                  outfmt=5) # must output XML!
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate() # forces the main script to wait

def local_tblastx_2file(query_file, dbfile_path, outfile, prefs):
    """Perform blastx against local database."""
    cline = NcbitblastxCommandline(query=query_file,
                                  db=dbfile_path,
                                  out=outfile,
                                  evalue=prefs['evalue'],
                                  outfmt=prefs['outfmt_pref'])
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate() # forces the main script to wait

def local_tblastn_2file(query_file, dbfile_path, outfile, prefs):
    """Perform blastx against local database."""
    cline = NcbitblastnCommandline(query=query_file,
                                  db=dbfile_path,
                                  out=outfile,
                                  evalue=prefs['evalue'],
                                  outfmt=prefs['outfmt_pref'])
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate() # forces the main script to wait

def local_blastx_2file(query_file, dbfile_path, outfile, prefs):
    """Perform blastx against local database."""
    cline = NcbiblastxCommandline(query=query_file,
                                  db=dbfile_path,
                                  out=outfile,
                                  evalue=prefs['evalue'],
                                  outfmt=prefs['outfmt_pref'])
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate() # forces the main script to wait

def local_rpsblast_2file(query_file, dbfile_path, outfile, prefs):
    """Perform RPS Blast against local database."""
    cline = NcbirpsblastCommandline(query=query_file,
                                  db=dbfile_path,
                                  out=outfile,
                                  evalue=prefs['evalue'],
                                  outfmt=5) # must output XML!
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate() # forces the main script to wait

def remote_blastp_2file(query_string, database, outfile, evalue):
    """Perform blastp against remote database."""
    result_handle = NCBIWWW.qblast('blastp',
                                   database,
                                   query_string,
                                   expect=evalue)
    save_file = open(outfile, "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()

def read_array(filename, dtype, separator='\t'):
    # From Numpy cookbook
    """ Read a file with an arbitrary number of columns.
        The type of data in each column is arbitrary
        It will be cast to the given dtype at runtime
    """
    cast = numpy.cast
    data = [[] for dummy in xrange(len(dtype))]
    for line in open(filename, 'r'):
        fields = line.strip().split(separator)
        for i, number in enumerate(fields):
            data[i].append(number)
    for i in xrange(len(dtype)):
        #print data[i]
        data[i] = cast[dtype[i]](data[i])
    return numpy.rec.array(data, dtype=dtype)

# Blast results arrays datatypes
blast_dtypes = numpy.dtype([('query', 'S16'),
                           ('dbhit', 'S32'),
                           ('idp', 'float'),
                           ('mlen', 'uint8'),
                           ('mms', 'uint8'),
                           ('gaps', 'uint8'),
                           ('q_start', 'uint32'),
                           ('q_end', 'uint32'),
                           ('r_start', 'uint32'),
                           ('r_end', 'uint32'),
                           ('evalue', 'S5'),
                           ('bitscore', 'float')])