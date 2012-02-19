# script to rename contigs in multifasta files

from genomes import all as genomes
from libs import load_multifasta, write_fasta

for genome in genomes:
    print genome['file']
    file_path = "data/genomes/"+genome['file']
    outfile_path = "data/renamed/"+genome['file']
    contigs = load_multifasta(file_path)
    renamed = []
    counter = 1
    for contig in contigs:
        contig.id = genome['name']+"_"+str(counter)
        contig_path = "data/contigs/"+contig.id+".fas"
        write_fasta(contig_path, contig)
        renamed.append(contig)
        counter +=1
    write_fasta(outfile_path, renamed)