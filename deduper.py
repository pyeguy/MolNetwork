from sys import argv
from collections import Counter
script, infile = argv

infile_components = infile.split('.')

# used_names = set()
counter = Counter()
deduped = []
with open(infile) as fin:
	for line in fin:
		line = line.strip()
		if line:
			line = line.split('\t')
			if line[1] not in counter:
				deduped.append(line)
			counter[line[1]] +=1

with open('{}_deduped.smi'.format(infile_components[0]),'w') as fout:
	for line in deduped:
		print('{}\t{}'.format(*line),file=fout)
with open('{}_counts.txt'.format(infile_components[0]),'w') as fout:
	for line in counter.items():
		print('{}\t{}'.format(*line),file=fout)
