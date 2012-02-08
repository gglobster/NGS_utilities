import re
from os import listdir
from Bio import SeqIO
from os import path, makedirs

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