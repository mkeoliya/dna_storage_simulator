import random
import argparse

DICT = 'AGTC'


def gen_sample_dj(gold, sub_p, del_p, ins_p):
    res = []
    for w in gold:
        r = random.random()
        if r < sub_p:
            res.append(random.choice(DICT))
        elif r < sub_p + ins_p:
            res.append(random.choice(DICT))
            res.append(w)
        elif r > sub_p+ins_p+del_p:
            res.append(w)
    return ''.join(res)

def simulate_strand(strand, cov, sub_p, del_p, ins_p):
    res = []
    for _ in range(cov):
        res.append(gen_sample_dj(strand, sub_p, del_p, ins_p))
    return res


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Noise generation')
    parser.add_argument('--N', type=int, default=5)
    parser.add_argument('--R', type=int, default=0)
    parser.add_argument('--subs', type=float, default=0.01)
    parser.add_argument('--dels', type=float, default=0.01)
    parser.add_argument('--inss', type=float, default=0.01)
    parser.add_argument("--i", type=str)
    parser.add_argument(
        "--o", type=str, default="synth_data_clustered.txt")

    args = parser.parse_args()
    f = open(args.i)
    strands = f.readlines()
    f.close()
    f = open(args.o, "w")
    i = 0
    for strand in strands:
        cluster = []
        f.write(strand.strip() + "\n")
        f.write('*****************************\n')
        i = i+1
        myN = args.N
        if args.R == 1:
            myN = random.randint(args.N/2, args.N*1.5)
        for x in range(0, myN):
            cluster.append(gen_sample_dj(
                strand.strip(), args.subs, args.dels, args.inss))
        for s in cluster:
            f.write(s.strip()+"\n")

        f.write("\n\n")
