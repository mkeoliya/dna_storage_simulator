
'''
Generate clustered data of the form:

<num>
noisy copy
noisy copy
noisy copy

<num>
noisy copy
noisy copy
'''
import random

with open('../../Data/microsoft-real/Centers.txt', 'r') as f:
    strands = f.readlines()

with open('../../Data/microsoft-real/Clusters.txt', 'r') as f2:
    clusters = f2.readlines()

cov5_f = open('../../Data/cov5/nanopore-real-noshuf/clusters.txt', "w")
cov6_f = open('../../Data/cov6/nanopore-real-noshuf/clusters.txt', "w")

cov5_refs = open('../../Data/cov5/nanopore-real-noshuf/refs.txt', "w")
cov6_refs = open('../../Data/cov6/nanopore-real-noshuf/refs.txt', "w")

i = 0
count = 0
c5, c6 = 0, 0
current_cluster = []
lengths = []
new_thing = []
for line in clusters:
    if '=' in line:
        if len(current_cluster) >= 6:
            new_thing.append(strands[i])
            cov5_refs.write(strands[i])
            cov6_refs.write(strands[i])
        
            cov5_f.write(str(5)+"\n")
            c5 += 1
            for s in current_cluster[:5]:
                cov5_f.write(s)
            
            cov6_f.write(str(6)+"\n")
            c6 += 1
            for s in current_cluster[:6]:
                cov6_f.write(s)

        current_cluster = []
        i += 1
    else:
        current_cluster.append(line)

print(len(new_thing))