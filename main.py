# Import Packages First
import math
import matplotlib.pyplot as plt
import itertools
import portion as P

if __name__ == '__main__':
    file_path = 'data/GB_GCA_000008085.1_protein.txt'


class EcodDomain:
    """Class for working with ECOD domains."""

    def __init__(self, ecod_line):
        self.uid = ecod_line.split('\t')[0]
        self.ecod_domain_id = ecod_line.split('\t')[1]
        self.f_id = ecod_line.split('\t')[3]
        self.pdb = ecod_line.split('\t')[4]
        self.chain = ecod_line.split('\t')[5]
        self.pdb_range = self.parse_ecod_range(ecod_line.split('\t')[6])
        self.seqid_range = ecod_line.split('\t')[7]
        self.arch_name = ecod_line.split('\t')[9]
        self.x_name = ecod_line.split('\t')[10]
        self.h_name = ecod_line.split('\t')[11]
        self.t_name = ecod_line.split('\t')[12]
        self.f_name = ecod_line.split('\t')[13]
        self.asm_status = ecod_line.split('\t')[14]
        self.ligand = ecod_line.split('\t')[15].replace('NO_LIGANDS_4A',
                                                        '').split(',')

    def __str__(self):
        """Prints string with some useful information"""
        return f"{self.ecod_domain_id}: {self.f_id} ({self.x_name}, {self.f_name})"

    def __len__(self):
        """Returns the length of the ECOD domain."""
        pass
        # TODO: create __len__

    def __eq__(self, other):
        """Compare to uid or ecod_domain_id with == operator."""
        if type(other) == EcodDomain:
            return self.uid == other.uid
        elif other == self.uid or other == self.ecod_domain_id:
            return True
        else:
            return False

    def __ne__(self, other):
        """Defines != operator; see __eq__."""
        return not self.__eq__(other)

    def __contains__(self, item):
        """check ecod hierarchy with in operator"""
        if item == self.f_id:
            return True
        elif self.f_id.startswith(item + '.'):
            return True
        else:
            return False

    def structure_path(self):
        """
        Return path to structure in ECOD structure F99 database.
        Does not check if structure Exists.
        """
        return f""
        # TODO: what is this?

    def parse_ecod_range(self, ecod_range):
        """Parse ECOD range."""
        pass
        # TODO: (8) create parse_ecod_range

    def pymol_selector(self):
        """Return PyMOL selection command for domain."""
        return f"select {self.ecod_domain_id}, {self.pdb} and chain {self.chain} and resi"


class HmmerHit:
    """Class for working with HMMER hits"""

    def __init__(self, hmm_line):
        hmm_line = hmm_line.split()

        self.hmm_name = hmm_line[3]
        self.hmm_length = int(hmm_line[5])
        self.hmm_coverage = (int(hmm_line[16]) - int(hmm_line[15]) + 1) / float(
            hmm_line[5])
        self.gene_name = hmm_line[0]
        self.gene_length = hmm_line[2]
        self.hmm_range = [int(hmm_line[15]), int(hmm_line[16])]
        self.hit_range = [int(hmm_line[19]),
                          int(hmm_line[20])]  # range of hit within gene
        self.evalue = float(hmm_line[6])
        self.ievalue = float(hmm_line[12])

    def __str__(self):
        """Prints string with some useful information"""
        return (f"{self.hmm_name}: len: {self.hmm_length}, "
                f"coverage: {self.hmm_coverage}, gene_name: {self.gene_name}, "
                f"gene_len: {self.gene_length}, hmm_range: {self.hmm_range}, "
                f"hit_range: {self.hit_range}, evalue: {self.evalue}, "
                f"ievalue: {self.ievalue}")

    def __eq__(self, other):
        """Compare the hmm_names of self and other using == operator."""
        if type(other) == HmmerHit:
            return self.hmm_name == other.hmm_name
        else:
            return False

    def __ne__(self, other):
        """Defines != operator; see __eq__."""
        return not self.__eq__(other)


def parseHmm(file_path):
    """
    Returns a dictionary of HmmerHits for each gene.
    Takes the file path to the HMMER output file as input.
    """

    with open(file_path, 'r') as file:
        output = file.readlines()

    gene_hits = {}

    for line in output:
        if line[0] == '#':  # skip the first lines
            continue
        hit = HmmerHit(line)
        if hit.gene_name not in gene_hits:
            gene_hits[hit.gene_name] = [hit]
        else:
            gene_hits[hit.gene_name].append(hit)
    return gene_hits


def overlap(range1, range2):
    '''
    Return overlap length between two ranges as integer.
    Return 0 if no overlap.
    Assumes [start, end] order for both range1 and 2.
    '''

    # Check that [start, end] order is correct
    if range1[0] > range1[1] or range2[0] > range2[1]:
        return None

    max_start = max(range1[0], range2[0])
    min_end = min(range1[1], range2[1])
    if max_start <= min_end:
        return min_end - max_start + 1
    else:
        return 0

def overlaP(range1, range2):
    '''
    Return True if two ranges overlap, False otherwise.
    Uses portion package.
    '''
    return P.closed(range1[0], range1[1]).overlaps(P.closed(range2[0], range2[1]))


# TODO: (5) finish reconciliation function
def fetch_representative_hits(hits_dict_list,
                              e_value_threshold,
                              length_threshold,
                              coverage_threshold):
    # First filtering ste
    filtered_hits = list(filter(lambda x: (x['ievalue'] < e_value_threshold
                                           and (x[
                                                    'hmm_coverage'] > coverage_threshold)),
                                hits_dict_list))

    # Second filtering step
    # .sort(key=lambda x: (x[0], -(x[-1] - x[-2])))
    # accepted_hits = [hits_list_copy[0]]
    # for hit in hits_list_copy[1:]:
    #     for ahit in accepted_hits:
    #         if overlap(ahit[1:], hit[
    #                              1:]) == 0:  # We will use something less stringent in final script
    #             accepted_hits.append(hit)




gene_hits = parseHmm(file_path)
# for gene, hits in gene_hits.items():
#     print(f'{gene}: {hits}')

for hit in gene_hits['AE017199.1_291']:
    print(hit.hmm_name)

for genes, hits in gene_hits.items():
    print(f'<<{genes}>>')
    for hit in hits:
        print(hit)
