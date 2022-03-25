import copy
import strand_error_sim
import os
import numpy as np

class Simulator:
    """
    Class variables:
    :var self.total_error_rates: Dictionary of the total error rates used in the simulation, as provided in error_rates
        parameter.
    :var self.base_error_rates: Error rates corresponding to each base, as passed.
    :var self.long_deletion_length_rates: Based on (excluding single base deletion):
        https://www.biorxiv.org/content/biorxiv/early/2019/11/13/840231/F12.large.jpg?width=800&height=600&carousel=1
    """
    def __init__(self, total_error_rates, base_error_rates, long_deletion_length_rates):
        """
        @param total_error_rates: Dictionary of the total error rates used in the simulation.
            Example of a dictionary:
            {'d': 0.1, 'i': 0.2, 's': 0.1, 'ld': 0.6}
            NOTE: Dictionary can be passed with values as strings, as the constructor can to parse them to floats.
        @param base_error_rates: Dictionary of dictionaries for each base.
            Example:
            {   'A': {'s': 0.1, 'i': 0.2, 'pi': 0.1, 'd': 0.05, 'ld': 0.6},
                'T': {...},
                'C': {...},
                'G': {...}
            }
        """
        self.total_error_rates = total_error_rates
        self.base_error_rates = base_error_rates
        self.long_deletion_length_rates = long_deletion_length_rates

    def simulate_strand(self, original_strand, num_copies):
        # for each strand, do the simulation on a copy of it num_copies (the generated number of copies) times:
        output_strands = []
        for j in range(num_copies):

            # duplicate strand to create a working (output) strand:
            output_strand = copy.deepcopy(original_strand)

            # create a strand simulator for it:
            strand_error_simulator = strand_error_sim.StrandErrorSimulation(self.total_error_rates,
                                                                            self.base_error_rates,
                                                                            self.long_deletion_length_rates,
                                                                            output_strand)
            # simulate 
            output_strand = strand_error_simulator.simulate_errors_on_strand()

            output_strands.append(output_strand)
        
        return output_strands
    
    def simulate_strands(self, original_strands, num_copies):
        return [self.simulate_strand(s, num_copies) for s in original_strands]

def get_pos_for_v_shape(skew, i, mid):
    return int((len(skew) - 1) * abs(i - mid) / mid)

def get_pos_for_a_shape(skew, i, mid):
    pos = -1 * int(len(skew) * abs(i - mid) / mid)
    if i == 55:
        pos = -1
    return pos

if __name__ == '__main__':
    total_error_rates = {
        'd': 0.0833,
        'i': 0.0833,
        's': 0.0833,
        'ld': 0,
    }
    
    base_error_rates = {
        'A': {
            'pi': 0.01588836102870499,
            'i': 0.01621475841752588,
            's': 0.015851018954560227,
            'd': 0.011764136395383413,
            'ld': 0.003110456472280425,
        },
        'G': {
            'pi': 0.023879556701749495,
            'i': 0.018091587091377988,
            's': 0.016569317627378545,
            'd': 0.017090680978601797,
            'ld': 0.003931135429009549
        },
        'C': {
            'pi': 0.017915924372554178,
            'i': 0.01887005911546257,
            's': 0.01853747500507736,
            'd': 0.01712126643667476,
            'ld': 0.0033067584089939467
        },
        'T': {
            'pi': 0.01581873288059774,
            'i': 0.020162401335515973,
            's': 0.017022155950663705,
            'd': 0.014674738110288119,
            'ld': 0.002907934758532013,
        }
    }

    long_deletion_length_rates = {2: 84,
                                  3: 13,
                                  4: 1.8,
                                  5: 0.2,
                                  6: 0.02}
    

    sim = Simulator(total_error_rates, base_error_rates, long_deletion_length_rates)

    skew = sorted(np.random.triangular(0, 0.15, 0.3, 20000))
    mid = 55

    skew_vals = [0] * 111
    for i in range(0, 111):
        pos = get_pos_for_a_shape(skew, i, mid)
        skew_vals[i] = skew[pos]

    file_path = os.path.dirname(__file__)
    if file_path != "":
        os.chdir(file_path)

    refs_file = open('../../Data/microsoft-real/Centers.txt', 'r')
    original_strands = [line.strip() for line in refs_file.readlines()]

    i = 0

    # cov_5 = open('../../Data/cov5/nanopore-2nd-order-skew/clusters.txt', 'w')
    # cov_6 = open('../../Data/cov6/nanopore-2nd-order-skew/clusters.txt', 'w')
    cov_10 = open('../../Data/cov10/profs-params/clusters.txt', 'w')
    cov_10_refs = open('../../Data/cov10/profs-params/refs.txt', 'w')


    strand_cov = 10
    example_strand = "TTGTCACTAGAGGACGCACGCTCTATTTTTATGATCCATTGATGTCCCTGACGCTGCAAAATTTGCAACCAGGCAGTCTTCGCGGTAGGTCCTAGCCTTATTGTCACTAGAGGACGCACGCTCTATTTTTATGATCCATTGATGTCCCTGACGCTGCAAAATTTGCAACCAGGCAGTCTTCGCGGTAGGTCCTAGCCTTA"
    for strand in original_strands:
        # cov_5.write(f'{5}\n')
        # cov_6.write(f'{6}\n')
        cov_10.write(f'{10}\n')

        new_str = (strand + strand)[:200]
        sims = sim.simulate_strand(new_str, strand_cov)
        # for s in sims[:5]:
        #     cov_5.write(s + '\n')
        
        # for s in sims[:6]:
        #     cov_6.write(s + '\n')
        cov_10_refs.write(new_str+'\n')
        for s in sims:
            cov_10.write(s + '\n')

        i += 1
