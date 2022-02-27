
'''
Generate clustered data of the form:

original 
*********
noisy copy
noisy copy
noisy copy

original 
**********
noisy copy
noisy copy
'''

with open('Centers.txt') as f:
    strands = f.readlines()

with open('Clusters.txt') as f2:
    clusters = f2.readlines()

i = 0
current_cluster = []
output = open('nanopore_clustered_filtered_n_10.txt', "w")
lengths = []
for line in clusters:
    if '=' in line:
        output.write(strands[i].strip() + "\n")
        output.write('*****************************\n')

        lengths.append(len(current_cluster))
        for s in current_cluster:
            output.write(s.strip()+"\n")
        output.write("\n\n")
        current_cluster = []
        i += 1
    else:
        current_cluster.append(line)

print(lengths)