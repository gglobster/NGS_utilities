## script to predict orfs on sequences

import re
from os import path
from sys import argv
from libs import from_dir, ensure_dir, fas2gbk, gbk2fas, write_genbank, \
    load_genbank, train_prodigal, run_prodigal, load_multifasta, \
    local_blastp_2file, collect_cogs
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import generic_dna

origin_dir = "data/"+argv[1]+"/"
seq_dir = origin_dir+argv[2]+"/"
file_ext = argv[3]

if len(argv) < 5:
    trim_ids = ''
else:
    trim_ids = argv[4]

blast_dir = origin_dir+"blast/"
prot_db = "data/ref_dbs/Bacteria_prot"

annot_gbk_dir = origin_dir+"annot_gbk/"
annot_aa_dir = origin_dir+"annot_aa/"
trn_file = origin_dir+"prodigal.trn"

ensure_dir([annot_gbk_dir, annot_aa_dir, blast_dir])

filenames = from_dir(seq_dir, re.compile(r'.*\.'+file_ext+'.*'))

for filename in filenames:
    rec_name = filename[:filename.find(trim_ids+"."+file_ext)]

    print rec_name, "...",

    # load data
    if file_ext == 'fas':
        fas_file = seq_dir+"/"+filename
        gbk_file = fas2gbk(fas_file)
        record = load_genbank(gbk_file)
    else:
        gbk_file = seq_dir+"/"+filename
        fas_file = gbk2fas(gbk_file)
        record = load_genbank(gbk_file)

    # run prediction
    annot_aa = annot_aa_dir+rec_name+"_ann.fas"
    annot_gbk = annot_gbk_dir+rec_name+"_ann.gbk"
    if not path.exists(trn_file):
        train_prodigal(fas_file, trn_file, "-q")
    if not path.exists(annot_aa):
        run_prodigal(fas_file, annot_gbk, annot_aa, trn_file, "-q")

    # blast the amino acids against COG
    blast_out = blast_dir+rec_name+".xml"
    blast_prefs = {'evalue': 0.01, 'outfmt_pref': 6}
    if not path.exists(blast_out):
        local_blastp_2file(annot_aa, prot_db, blast_out, blast_prefs)
    # collect best hits
    rec_cogs = collect_cogs(blast_out)

    # collect orfs
    record.features = []
    aa_record = load_multifasta(annot_aa)
    counter = 1
    for aa_rec in aa_record:
        this_prot = 'Query_'+str(counter)
        annotation = rec_cogs[this_prot]
        # get feature details from description line
        # because prodigal output fails to load as valid genbank
        defline = aa_rec.description
        pattern = re.compile('.+#\s(\d+)\s#\s(\d+)\s#\s(\S*1)\s#\sID.+')
        match = pattern.match(defline)
        start_pos = int(match.group(1))
        end_pos = int(match.group(2))
        strand_pos = int(match.group(3))
        feat_loc = FeatureLocation(start_pos-1, end_pos) # adjust for 0-index
        l_tag = rec_name+"_"+str(counter)
        # consolidation feature annotations
        quals = {'note': defline, 'locus_tag': l_tag,
                 'fct': annotation, 'translation': aa_rec.seq}
        feature = SeqFeature(location=feat_loc,
                             strand=strand_pos,
                             id='cds_'+str(counter),
                             type='CDS',
                             qualifiers=quals)
        record.features.append(feature)
        counter +=1

    # add annotations for Nx100 spacers
    sequence = str(record.seq)
    separator = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
    regex = re.compile(separator, re.IGNORECASE)
    spacers = [match.start() for match in regex.finditer(str(record.seq))]
    #print spacers
    for spacer in spacers:
        space_loc = FeatureLocation(spacer, spacer+100)
        feature = SeqFeature(location=space_loc,
                             type='spacer')
        record.features.append(feature)

    # save record with annotations
    record.description = rec_name+"_with_ORFs"
    record.name = rec_name
    record.dbxrefs = ["Project: "+argv[1]+"/"+rec_name]
    record.seq.alphabet = generic_dna
    write_genbank(annot_gbk, record)

    print "OK"


