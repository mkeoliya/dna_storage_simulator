import os
import naive_sim

def set_file_path():
    file_path = os.path.dirname(__file__)
    if file_path != "":
        os.chdir(file_path)

if __name__ == '__main__':
    
    total_error_rates = {
        'd': 0.015163160437436866,
        'i': 0.018348000570951734,
        's': 0.01700159756686723 ,
    }

    set_file_path()
    refs_file = open('../../Data/microsoft-real/Centers.txt', 'r')
    original_strands = [line.strip() for line in refs_file.readlines()]

    i = 0
    cov_5 = open('../../Data/cov5/nanopore-sim/clusters.txt', 'w')
    cov_6 = open('../../Data/cov6/nanopore-sim/clusters.txt', 'w')

    for strand in original_strands:
        cov_5.write(f'{5}\n')
        cov_6.write(f'{6}\n')

        coverage = 6
        sims = naive_sim.simulate_strand(strand, coverage, total_error_rates['s'], total_error_rates['d'], total_error_rates['i'])

        for s in sims[:5]:
            cov_5.write(s + '\n')
        
        for s in sims[:6]:
            cov_6.write(s + '\n')

        i += 1
