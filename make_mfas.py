## script to combine several fasta files into a single one

import re
from sys import argv
from libs import from_dir, load_fasta, write_fasta

origin_dir = "data/"+argv[1]
destin_file = origin_dir+"/"+argv[2]+".fas"

filenames = from_dir(origin_dir, re.compile(r'.*\.fas.*'))

records = []

for filename in filenames:
	# load record 
	records.append(load_fasta(origin_dir+"/"+filename))
	
	print filename

write_fasta(destin_file, records)