# script to do batch blast

from os import path
from sys import argv
from libs import make_blastDB, local_tblastn_2file, load_genbank, write_fasta

from genomes import all as genome_list

data_dir = "data/"+argv[1]+"/"
infile = data_dir+argv[2]

# for record in list:
for genome in genome_list:
    print genome['name'],
    # mame a blast DB if it doesn't already exist
    genome_path = "data/genomes/"+genome['file']
    dbfile_path = "data/blast_db/"+genome['name']
    if genome['input'] == 'cgbk':
        print "ignoring cgbk file"
    else:
        if genome['input'] == 'gbk':
            try:
                record = load_genbank(genome_path)
                genome_path = "data/genomes/"+genome['name']+".fas"
                write_fasta(genome_path, record)
            except IOError:
                print "no file"
        if not path.exists(dbfile_path+".nhr"):
            try:
                make_blastDB(dbfile_path, genome_path, 'nucl')
            except IOError:
                print "failed to make DB"
        try:
            # blastx against each genome DB
            outfile = data_dir+genome['name']+".txt"
            prefs = {'evalue': 0.001, 'outfmt_pref': 6}
            local_tblastn_2file(infile, dbfile_path, outfile, prefs)
        except Exception:
            print "failed to blast"
        else:
            print "OK"
    