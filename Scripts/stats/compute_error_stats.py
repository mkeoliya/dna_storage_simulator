import random
import matching_new


class Base:
    def __init__(self, char):
        self.char = char
        self.total = 0
        self.subs = 0
        self.pre_ins = 0  # some base inserted before this character
        self.ins = 0  # this character inserted before some base
        self.dels = 0
        self.long_dels = 0

    def __repr__(self):
        return f'Total: {self.total}, Pre-Insertions: {self.pre_ins}, Insertions: {self.ins}, Subs: {self.subs}, Dels: {self.dels}, Long Dels: {self.long_dels}'

    def error_rates(self):
        return f'{self.char},\nPre-Insertions: {self.pre_ins / self.total},\nInsertions: {self.ins / self.total},\nSubs: {self.subs / self.total},\nDels: {self.dels / self.total},\nLong Dels: {self.long_dels / self.total}'


def read_strands():
    with open('../data/microsoft-real/Centers.txt') as f:
        return f.readlines()

def read_clusters():
    with open('../data/microsoft-real/Clusters.txt') as f:
        return f.readlines()


def count_ops(strand, read, ops, bases):
    strand_pos = 0
    read_pos = 0
    for op in ops:
        if strand_pos >= len(strand):  # trailing insertions are ignored
            continue

        base = strand[strand_pos]

        if op == '+':
            bases[base].dels += 1
        elif op == 's':
            bases[base].subs += 1
        elif op == '-':
            copy_base = read[read_pos]
            bases[base].pre_ins += 1
            bases[copy_base].ins += 1

        if op != '-':
            strand_pos += 1

        if op != '+':
            read_pos += 1


def count_bases(read, bases):
    for base in read:
        bases[base].total += 1


def count_long_dels(strand, ops, bases):
    long_deletion_length_rates = {2: 2.8 * (10 ** (-4)),
                                  3: 7.75 * (10 ** (-5)),
                                  4: 3.25 * (10 ** (-5)),
                                  5: 10 ** (-6),
                                  6: 5.5 * (10 ** (-8))}

    consec_dels = []
    strand_pos = 0
    for op in ops:
        if op == '+':
            consec_dels.append(strand[strand_pos])
        else:  # no more consecutive deletes
            if len(consec_dels) >= 2:
                actual_dels.append(len(consec_dels))
                while len(consec_dels) >= 2:
                    filtered_dict = {
                        k: v for k, v in long_deletion_length_rates.items() if k <= len(consec_dels)}
                    options = list(filtered_dict.keys())
                    rates = list(filtered_dict.values())
                    draw = random.choices(options, weights=rates, k=1)

                    deletion_length = draw[0]
                    bases[consec_dels[0]].long_dels += 1
                    for i in range(deletion_length):
                        if i < len(consec_dels):
                            bases[consec_dels[i]].dels -= 1
                    consec_dels = consec_dels[deletion_length:]

            consec_dels = []

        if op != '-':
            strand_pos += 1


def compute_error(data):
    bases = {
        'A': Base('A'),
        'G': Base('G'),
        'C': Base('C'),
        'T': Base('T')
    }

    for strand, cluster in data:
        for read in cluster:
            _, ops = matching_new.ops_list(strand, read)
            count_bases(read, bases)
            count_ops(strand, read, ops, bases)
            count_long_dels(strand, ops, bases)

    total = sum([b.total for b in bases.values()])
    e_ins = sum([b.pre_ins for b in bases.values()]) / total
    e_subs = sum([b.subs for b in bases.values()]) / total
    e_dels = sum([b.dels for b in bases.values()]) / total
    e_long_dels = sum([b.long_dels for b in bases.values()]) / total

    return e_ins, e_subs, e_dels, e_long_dels, bases


def read_data():
    strands = read_strands()
    clusters = read_clusters()

    res = [(s.strip(), []) for s in strands]
    i = 0
    for line in clusters:
        if '=' in line:
            i += 1
        else:
            _, cluster = res[i]
            cluster.append(line.strip())

    return res

actual_dels = []
if __name__ == '__main__':
    data = read_data()  # list of (ref, [reads])

    f = open('matching_ops.txt', 'w')
    for strand, cluster in data:
        for read in cluster:
            _, ops = matching_new.ops_list(strand, read)
            f.write(''.join(ops) + '\n')


    # e_ins, e_subs, e_dels, e_long_dels, bases = compute_error(data[:5])
    # print(e_ins, e_subs, e_dels, e_long_dels)
    # print(bases['A'].error_rates())
    # print(bases['G'].error_rates())
    # print(bases['C'].error_rates())
    # print(bases['T'].error_rates())
