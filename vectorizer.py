# script to represent sequences as vectors (based on gene annots)

from sys import argv
from libs import load_genbank, write_fasta, make_blastDB, \
    local_tblastn_2file, read_array, blast_dtypes, ensure_dir
from Bio.SeqRecord import SeqRecord

from sets.phages_set import all as genomes

data_dir = 'data/'+argv[1]+'/'
seq_dir = data_dir+argv[2]+'/'
out_dir = data_dir+argv[3]+'/'
feat_type = argv[4]
threshold = argv[5]

ensure_dir([out_dir])

db_file = out_dir+'ref_DB.fas'
db_path = out_dir+'refs'

new_DB = True
init_DB = False

symbolDB = {}
vectorDB = []

sym_cnt = 0

for genome in genomes:

    g_vector = []

    print genome['name'],

    while True:

        try:
            assert genome['input'] == 'gbk'
        except ValueError:
            print "bad format (skipping)"
            break

        # load genome file to extract features (to proteins in mfas file)
        record = load_genbank(seq_dir+genome['file'])
        select = [feat for feat in record.features if feat.type == feat_type]
        feat_cnt = 0

        # cycle through selected features
        for feat in select:

            feat_cnt +=1
            rec = feat.extract(record)
            rec.description = genome['name']+'_'+feat_type+'_'+str(feat_cnt)

            # initialize or update blast DB
            if init_DB:
                ref_records = [value[0] for value in symbolDB.values()]
                write_fasta(db_file, ref_records)
                try:
                    make_blastDB(db_path, db_file, 'nucl')
                except Exception:
                    print "failed to make blast DB"
                    exit()
                init_DB = False

            # first go: add all features as new symbols
            if new_DB:
                sym_cnt +=1
                symbol = 'S'+str(sym_cnt)
                rec.id = symbol
                symbolDB[symbol] = [rec]
                g_vector.append(symbol)

            else:
                # tblastn against the reference DB
                infile = data_dir+'temp.fas'
                outfile = data_dir+'temp.txt'
                write_fasta(infile, SeqRecord(rec.seq.translate(), id='temp'))
                prefs = {'evalue': 0.001, 'outfmt_pref': 6}
                try:
                    local_tblastn_2file(infile, db_path, outfile, prefs)
                except Exception:
                    print "failed to blast"
                    exit()

                # parse output -- take only first hit
                try:
                    hit = read_array(outfile, blast_dtypes)[0]
                except IndexError:
                    print '-',
                    # add record to DB as new symbol
                    sym_cnt +=1
                    symbol = 'S'+str(sym_cnt)
                    rec.id = symbol
                    symbolDB[symbol] = [rec]
                    g_vector.append(symbol)
                    init_DB = True
                else:
                    if hit[2] > int(threshold):
                        print '+',
                        # add record to DB as secondary match
                        symbol = hit[1]
                        rec.id = symbol+'_'+str(len(symbolDB[symbol])+1)
                        symbolDB[symbol].append(rec)
                        # add symbol to genome vector
                        g_vector.append(symbol)
                        init_DB = False
                    else:
                        print '<',
                        # add record to DB as new symbol
                        sym_cnt +=1
                        symbol = 'S'+str(sym_cnt)
                        rec.id = symbol
                        symbolDB[symbol] = [rec]
                        g_vector.append(symbol)
                        init_DB = True

        break

    vectorDB.append(g_vector)

    if new_DB:
        init_DB = True
        new_DB = False

    print 'OK'

# do something with the vectors
core_genome = [value for value in symbolDB.values()
               if len(value) == len(genomes)]
core_keys = [key for key in symbolDB
             if len(symbolDB[key]) == len(genomes)]

# output protein sequences by feature
for nt_set in core_genome:
    aa_set = [SeqRecord(nt_rec.seq, id=nt_rec.id+'_aa',
                        description=nt_rec.description) for nt_rec in nt_set]
    write_fasta(out_dir+nt_set[0].id+'_aa.fas', aa_set)

# concatenate aa seqs for later phylo analysis
counter = 0
concat_recs = []
while counter < len(genomes):
    nt_recs = [sym_set[counter] for sym_set in core_genome]
    aa_recs = [SeqRecord(nt_rec.seq.translate(), id=nt_rec.id+'_aa',
                         description=nt_rec.description) for nt_rec in nt_recs]

    aa_concat = SeqRecord('', '')
    for aa_rec in aa_recs:
        aa_concat += aa_rec

    aa_concat.id = genomes[counter]['name']
    aa_concat.description = ':'.join(core_keys)

    concat_recs.append(aa_concat)
    counter +=1

write_fasta(out_dir+'concat_aa.fas', concat_recs)

# TODO: must fix annotation coordinates!