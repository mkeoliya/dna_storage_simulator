
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
import os

file_path = os.path.dirname(__file__)
if file_path != "":
    os.chdir(file_path)

with open('../../data/microsoft-simulated/skewed-distribution/refs_part.txt') as f:
    strands = f.readlines()

with open('../../data/microsoft-simulated/skewed-distribution/recons-iter-part.txt') as f2:
    clusters = f2.readlines()


f = open('iter-skewed-simulated.csv', "w")
writer = csv.writer(f)
header = ['Refs', 'Reads']
writer.writerow(header)

for ref, cluster in zip(strands, clusters):
    row = [ref.strip(), cluster.strip()]
    writer.writerow(row)

