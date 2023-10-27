# Import Packages First
import math
import matplotlib.pyplot as plt
import itertools
import portion as P


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

    def parse_ecod_range(self, ecod_range):
        """Parse ECOD range."""
        pass

    def pymol_selector(self):
        """Return PyMOL selection command for domain."""
        return f"select {self.ecod_domain_id}, {self.pdb} and chain {self.chain} and resi"


class HmmerHit:
    """Class for working with HMMER hits"""
    pass


def parseHmm(file_path):
    # create gene_hits dictionary
    gene_hits = {}

    with open(file_path, 'r') as file:
        output = file.readlines()

    for line in output:
        if line[0] == '#':  # skip the first lines
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
    hmm_dict['hmm_coverage'] = (int(line[16]) - int(line[15]) + 1) / float(
        line[5])
    hmm_dict['gene_name'] = line[0]
    hmm_dict['gene_length'] = line[2]
    hmm_dict['hmm_range'] = [int(line[15]), int(line[16])]
    hmm_dict['gene_range'] = [int(line[19]),
                              int(line[20])]  # hit range, not glen
    hmm_dict['evalue'] = float(line[6])
    hmm_dict['ievalue'] = float(line[12])
    return hmm_dict


def overlapLength(range1, range2):
    """
    Return overlap length between two ranges as integer.
    Return 0 if no overlap.
    Return None if ranges are not in [start, end] order.

    Uses portion package.
    P.closed(0, 2) & P.closed(1, 3)
    [1,2]
    P.closed(0, 2) & P.closed(2, 3)
    [2]
    P.closed(0, 2) & P.closed(3, 4)
    ()
    """
    # Check that [start, end] order is correct
    if range1[0] > range1[1] or range2[0] > range2[1]:
        return None

    # return 0 if no overlap
    if not overlapBool(range1, range2):
        return 0

    overlap_range = overlapRange(range1, range2)
    if overlap_range[0] == overlap_range[1]:  # overlap by one
        return 1
    else:
        return overlap_range[1] - overlap_range[0] + 1


def overlapBool(range1, range2):
    """
    Return True if two ranges overlap, False otherwise.
    Return None if ranges are not in [start, end] order.

    Uses portion package.
    P.closed(1, 2).overlaps(P.closed(2, 3))
    True
    P.closed(1, 2).overlaps(P.open(2, 3))
    False
    """
    # Check that [start, end] order is correct
    if range1[0] > range1[1] or range2[0] > range2[1]:
        return None

    return P.closed(range1[0], range1[1]).overlaps(
        P.closed(range2[0], range2[1]))


def overlapRange(range1, range2):
    """
    Return a list [overlap_start, overlap_end] if two ranges overlap
    Return an empty list [] if no overlap.
    Return None if ranges are not in [start, end] order.

    Uses portion package:
    P.closed(0, 2) & P.closed(1, 3)
    [1,2]
    P.closed(0, 2) & P.closed(2, 3)
    [2]
    P.closed(0, 2) & P.closed(3, 4)
    ()
    """
    # Check that [start, end] order is correct
    if range1[0] > range1[1] or range2[0] > range2[1]:
        return None

    # return [] if no overlap
    if not overlapBool(range1, range2):
        return []
    else:
        overlap_range = P.closed(range1[0], range1[1]) & P.closed(range2[0], range2[1])
        return [overlap_range.lower, overlap_range.upper]



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


if __name__ == '__main__':
    file_path_ = 'GB_GCA_000008085.1_protein.txt'
