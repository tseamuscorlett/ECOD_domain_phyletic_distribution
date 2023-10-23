# Import Packages First
import math
import matplotlib.pyplot as plt
import itertools


class EcodDomain:
    """Class for working with ECOD Domains"""

    def __init__(self, ecod_line):
        self.uid = ecod_line[0]
        self.ecod_domain_id = ecod_line[1]
        self.f_id = ecod_line[3]
        self.x_id = ecod_line[3].split('.')[0]
        self.pdb = ecod_line[4]
        self.chain = ecod_line[5]
        self.pdb_range = self.parse_ecod_range(ecod_line[6])
        self.seqid_range = ecod_line[7]
        self.arch_name = ecod_line[9]
        self.x_name = ecod_line[10]
        self.h_name = ecod_line[11]
        self.t_name = ecod_line[12]
        self.f_name = ecod_line[13]
        self.asm_status = ecod_line[14]
        self.ligand = ecod_line[15].replace('NO_LIGANDS_4A', '').split(',')

    def __str__(self):
        return f"{self.ecod_domain_id}: {self.f_id} ({self.x_name})"

    def __len__(self):
        pass

    def structure_path(self):
        pass

    def parse_ecod_range(self, ecod_range):
        pass


def parse_hmm_output(file_path):
    # create gene_hits
    gene_hits = {}
    # trolololololo
    with open(file_path, 'r') as file:
        output = file.readlines()

    for line in output:
        if line[0] == '#':
            continue
        hit = create_hmm_dictionary(line)
        if hit['gene_name'] not in gene_hits:
            gene_hits[hit['gene_name']] = [hit]
        else:
            gene_hits[hit['gene_name']].append(hit)


def create_hmm_dictionary(line):
    hmm_dict = {}
    line = line.split()
    hmm_dict['hmm_name'] = line[3]
    hmm_dict['hmm_length'] = int(line[5])
    hmm_dict['hmm_coverage'] = (int(line[16]) - int(line[15]) + 1) / float(line[5])
    hmm_dict['gene_name'] = line[0]
    hmm_dict['gene_length'] = line[2]
    hmm_dict['hmm_range'] = [int(line[15]), int(line[16])]
    hmm_dict['gene_range'] = [int(line[19]), int(line[20])]  # hit range, not glen
    hmm_dict['evalue'] = float(line[6])
    hmm_dict['ievalue'] = float(line[12])
    return hmm_dict


def overlap(range1, range2):
    '''
    Return overlap between two ranges as integer.
    Assumes [start, end] order for both range1 and 2.
    '''

    # Check that data satisfies expectation
    if range1[0] > range1[1] or range2[0] > range2[1]:
        return None

    max_start = max(range1[0], range2[0])
    min_end = min(range1[1], range2[1])
    if max_start <= min_end:
        return min_end - max_start + 1
    else:
        return 0


def fetch_representative_hits(hits_dict_list,
                              e_value_threshold,
                              length_threshold,
                              coverage_threshold):
    # First filtering ste
    filtered_hits = list(filter(lambda x: (x['ievalue'] < e_value_threshold
                                           and (x['hmm_coverage'] > coverage_threshold)),
                                hits_dict_list))

    # Second filtering step
    # .sort(key=lambda x: (x[0], -(x[-1] - x[-2])))
    # accepted_hits = [hits_list_copy[0]]
    # for hit in hits_list_copy[1:]:
    #     for ahit in accepted_hits:
    #         if overlap(ahit[1:], hit[
    #                              1:]) == 0:  # We will use something less stringent in final script
    #             accepted_hits.append(hit)


if __name__ == '__main__':
    file_path_ = 'GB_GCA_000008085.1_protein.txt'
