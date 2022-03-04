
'''
Generate data of the form:

ref1, read1
ref1, read2
...
ref2, read1
ref2, read2
...
'''

import csv

with open('refs_microsoft.txt') as f:
    refs = f.readlines()

with open('recons-our-bma.txt') as f2:
    reads = f2.readlines()


f = open('our-bma.csv', "w")
writer = csv.writer(f)
header = ['Refs', 'Reads']
writer.writerow(header)

for (ref, read) in zip(refs, reads):
    row = [ref.strip(), read.strip()]
    writer.writerow(row)
