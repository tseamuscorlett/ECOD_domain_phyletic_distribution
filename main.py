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
        self.hmm_coverage = (int(hmm_line[16]) - int(hmm_line[15]) + 1) / float(hmm_line[5]) # Cannot be a range of 0
        self.gene_name = hmm_line[0]
        self.gene_length = hmm_line[2]
        self.hmm_range = P.closed(int(hmm_line[15]), int(hmm_line[16]))  # Cannot be len 0, portion object Interval
        self.hit_range = P.closed(int(hmm_line[19]), int(hmm_line[20]))  # Envelope coordinates; portion object Inteval
        self.evalue = float(hmm_line[6])
        self.ievalue = float(hmm_line[12])
        self.cevalue = float(hmm_line[11]) # double check

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
            raise None # be careful about combining cross-class and != within class results

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


def overlapLength(po1, po2):
    """
    Takes two portion Intervals
    Returns the length of the overlap as integer
    """
    z = po1 & po2
    if z.empty:
        return 0
    else:
        return z.upper - z.lower + 1


def overlapBool(po1, po2):
    """
    Takes two portion Intervals
    Returns True if they overlap, False if they don't
    A function for when you don't want to use the built-in ".overlaps" operator
    """
    return po1.overlaps(po2)


def overlapRange(po1, po2):
    """
    Takes two portion Intervals
    Returns the overlap as a portion Interval
    A function for when you don't want to use the built-in "&" operator
    """
    return po1 & po2


def overlapPercentage(po1, po2):
    """
    Takes two portion Intervals
    Returns the percentages for the two overlap as [float, float]
    """
    z = po1 & po2
    if z.empty:
        return [0, 0]
    else:
        overlap_len = z.upper - z.lower + 1
        return [float(overlap_len/(po1.upper - po1.lower + 1)), float(overlap_len/(po2.upper - po2.lower + 1))]


def fetch_representative_hits(hmm_hits_list,
                              e_value_threshold,
                              coverage_threshold):
    """
    Takes a list of hmm_hits of class HmmerHit
    Returns a list of accepted HMM hits above the custom e_value_threshold,
    coverage_threshold, and no overlap more than 25% between each hit
    Returns None if all hits get rejected by the thresholds
    """

    # Filtering step
    filtered_hits = list(filter(lambda x: (x.ievalue < e_value_threshold
                                           and (x.hmm_coverage > coverage_threshold)),
                                hmm_hits_list))

    # Sort by i-evalue
    filtered_hits.sort(key=lambda x: x.ievalue)
    if len(filtered_hits) == 0:
        return None
    else:
        accepted_hits = [filtered_hits[0]]

    # Reconciliation step
    for hit in filtered_hits[1:]:
        overlap_c = 0
        for ahit in accepted_hits:
            if overlapLength(ahit.hit_range, hit.hit_range) != 0:
                percentages = overlapPercentage(ahit.hit_range, hit.hit_range)
                if percentages[0] >= 0.25 or percentages[1] >= 0.25:
                    overlap_c += 1
        if overlap_c == 0:
            accepted_hits.append(hit)

    return accepted_hits


"""
# Testing print statements
"""

# gene_hits = parseHmm(file_path)
# domain_lengths = []
# hmm_lengths = []
# for gene, hits in gene_hits.items():
#     print(f'{gene}: {hits}')
#
# for hit in gene_hits['AE017199.1_291']:
#     print(hit.hmm_name)
#
# for genes, hits in gene_hits.items():
#     print(f'<<{genes}>>')
#     for hit in hits:
#         print(hit)
#         domain_lengths.append(hit.hit_range.upper - hit.hit_range.lower + 1)
#         hmm_lengths.append(hit.hmm_length)
# print(sum(domain_lengths)/len(domain_lengths))  # 103.3272921108742
# print(sum(hmm_lengths)/len(hmm_lengths))  # 163.89339019189765
# print(len(gene_hits))  # 382

"""
# Testing overlap reconciliation algorithm for one gene
"""
# gene_hits = parseHmm(file_path)
# hits_to_reconcile = gene_hits['AE017199.1_291']
#
# reconcile = []
# for hit in hits_to_reconcile:
#     reconcile.append(hit.hmm_name)
# print(f'To reconcile: {reconcile}')
#
# accepted = fetch_representative_hits(hits_to_reconcile, 1e-5, 0.6)
# for hit in accepted:
#     print(hit.hmm_name)

"""
# Applying overlap reconciliation algorithm for WHOLE GENOME
"""


def wholeGenomeOverlapRecon(gene_hits, e_value_threshold,
                              coverage_threshold):
    """
    Returns a dictionary {gene: [HmmerHit]} after filtering &
    overlap reconciliation for the whole genome "gene_hits", a parseHmm output {gene: [HmmerHit]}

    Removes genes for which all the hits have been removed {gene: None}
    """
    result = {}
    for gene, hits in gene_hits.items():
        hits_to_keep = fetch_representative_hits(hits, e_value_threshold, coverage_threshold)
        if hits_to_keep is not None:
            result[gene] = hits_to_keep
    return result


# given an HMM_name, can you search if this genome has that HMM_profile?
def domainSearch(domain_name, accepted_hits):
    """
    Returns True if domain_name exists within the list of accepted_hits
    Returns False otherwise
    """
    for gene, hits in accepted_hits.items():
        for hit in hits:
            if hit.hmm_name == domain_name:
                return True
    return False

# parse HMMER output
gene_hits = parseHmm(file_path)

# reconcile overlaps
result = wholeGenomeOverlapRecon(gene_hits, 1e-5, 0.6)
print(result)

# search for domains
print(domainSearch('MCM_AAA', result))
print(domainSearch('fake_domain_name', result))


