import re
from os import listdir
from Bio import SeqIO

def load_multifasta(seqfile):
    """Load multiple records from Fasta file."""
    input_handle = open(seqfile, 'rU')
    multifasta = SeqIO.parse(input_handle, 'fasta')
    fasta_list = list(multifasta)
    for record in fasta_list:
        assert record.id
    return fasta_list

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