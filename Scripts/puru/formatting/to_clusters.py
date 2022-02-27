with open('SimulatedClusters.txt') as f:
    lines = f.readlines()

with open('Centers-9999.txt') as f:
    refs = f.readlines()

refs_file = open('refs_microsoft_simulated.txt', 'w')
cluster_file = open('clusters_microsoft_simulated.txt', 'w')


def write_cluster(i, ref, cluster):
    cluster_file.write(f'{len(cluster)}\n')
    for read in cluster:
        cluster_file.write(read)
    refs_file.write(ref)

ref_pos = 0
cluster_pos = 0
curr = []
for line in lines:
    if '=' in line:
        if len(curr) > 0:
            write_cluster(cluster_pos, refs[ref_pos], curr)
            cluster_pos += 1
        ref_pos += 1
        curr = []
    else:
        curr.append(line)

# last cluster
write_cluster(cluster_pos, refs[ref_pos], curr)

refs_file.close()
cluster_file.close()
