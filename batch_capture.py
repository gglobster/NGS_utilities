# script to capture sequences from the results of a batch blast

from sys import argv
from Bio.SeqRecord import SeqRecord
from libs import read_array, blast_dtypes, load_fasta, write_fasta

data_dir = "data/"+argv[1]+"/"
main_in = data_dir+argv[2]+"_results.txt"
main_out = data_dir+argv[2]+"_ctxt.fas"
capture_span = int(argv[3])

records = []

rec_array = read_array(main_in, blast_dtypes)

for line in rec_array:

    query = line[0]
    subject = line[1]

    rev_flag = False
    if line[8] < line[9]:
        q_start, q_stop = line[8]-1, line[9]
        rev_flag = False
    else:
        q_start, q_stop = line[9]-1, line[8]
        rev_flag = True

    c_start, c_stop = q_start-capture_span, q_stop+capture_span

    master_seq = load_fasta("data/contigs/"+subject+".fas")

    if c_start < 0:
        c_start = 0
    if c_stop > len(master_seq.seq):
        c_stop = len(master_seq.seq)

    seq_bit = master_seq[c_start:c_stop]

    if rev_flag:
        seq_bit = seq_bit.reverse_complement()
    record = SeqRecord(id=subject, seq=seq_bit.seq)
    records.append(record)

    rec_file = data_dir+subject+"_"+query+"_ctxt.fas"
    write_fasta(rec_file, record)

write_fasta(main_out, records)
