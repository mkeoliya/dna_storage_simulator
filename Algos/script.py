with open('EncodedStrands.txt') as f:
    refs = f.readlines()

with open('datap012/datap0.12.txt') as f:
    reads = f.readlines()


output = open('output.txt', 'w')
i = 0
for ref in refs:
    output.write(ref)
    output.write('*****************************\n')
    i += 1 # line containing number
    for _ in range(10):
        output.write(reads[i])
        i += 1
    output.write('\n\n')



