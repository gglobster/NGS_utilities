## script to unpack a single fasta file into separate ones

import re
from sys import argv
from libs import from_dir, load_multifasta, write_fasta

origin_dir = "data/"+argv[1]+"/"
origin_file = origin_dir+argv[2]+".fas"
destin_dir = "data/"+argv[3]+"/"

records = load_multifasta(origin_file)

for record in records:
	write_fasta(destin_dir+record.id+".fas", record)
	print record.id
	
