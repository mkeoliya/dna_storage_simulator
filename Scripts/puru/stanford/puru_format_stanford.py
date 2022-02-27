with open('file.txt') as f:
    lines = f.readlines()


refs_file = open('refs_stanford.txt', 'w')
clust_file = open('clusters_stanford.txt', 'w')

i = 0
clust_num = 0
while True:
    if i >= len(lines):
        break

    num_reads = int(lines[i].strip())

    if num_reads >= 10:
        i += 1
        # first line in group is ref
        refs_file.write(lines[i])

        # now, write reads
        clust_file.write(f'CLUSTER {clust_num}\n')
        for j in range(num_reads):
            i += 1
            if j >= 10:
                continue
            clust_file.write(lines[i])

        clust_num += 1
    else:
        i += num_reads + 1

    i += 1

refs_file.close()
clust_file.close()
