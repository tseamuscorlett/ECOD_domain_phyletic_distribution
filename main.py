# Import Packages First
import csv
import math
import os

import matplotlib.pyplot as plt
import itertools
import pandas as pd
import re
import portion as P

if __name__ == '__main__':
    file_path = 'data/GB_GCA_000008085.1_protein.faa'


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
        self.cevalue = float(hmm_line[11])  # double check

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
            raise None  # be careful about combining cross-class and != within class results

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
        return [float(overlap_len / (po1.upper - po1.lower + 1)),
                float(overlap_len / (po2.upper - po2.lower + 1))]


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
    Returns the frequency of domain_name within the list of accepted_hits as int
    """
    count = 0
    for gene, hits in accepted_hits.items():
        for hit in hits:
            if hit.hmm_name == domain_name:
                count += 1
    return count


"""
# Testing domainSearch
"""
# # parse HMMER output
# gene_hits = parseHmm(file_path)
#
# # reconcile overlaps
# result = wholeGenomeOverlapRecon(gene_hits, 1e-5, 0.6)
# print(result)
# print(len(result))  # 298 genes
#
# # print frequency of domains
# print(domainSearch('MCM_AAA', result))
# print(domainSearch('fake_domain_name', result))
#
# counter = 0
# for gene, hits in result.items():
#     counter += len(hits)
# print(counter)  # 414 hits; some domains may exist in multiple genes


"""
# Let's make a dataframe
"""


def populateDataFrame(file_path, df):
    """
    Parse the input faa file of Hmmer result, turns it into a "reconciled"
    dictionary of {gene: [HmmerHit]} using wholeGenomeOverlapRecon.
    Creates one "row" of dataframe where the row_name = genome name and
    col_names = domain names.
    Returns a new dataframe where the input df & the new "row" is concatenated.
    """

    # create a dict {gene: [HmmerHit]} from file_path
    gene_hits = parseHmm(file_path)
    # extract genome name
    pattern = '[A-Z]{2}_[A-Z]{3}_\d+\.\d+'
    genome_name = re.search(pattern, file_path).group()

    # create a "reconciled" dict {gene: [HmmerHit]}
    reconciled = wholeGenomeOverlapRecon(gene_hits, 1e-5, 0.6)

    # Traverse reconciled hits; accumulate domainSearch results
    searched = []
    new_data = {}
    for gene, hits in reconciled.items():
        for hit in hits:
            # new domain
            if hit.hmm_name not in searched:
                new_data[hit.hmm_name] = domainSearch(hit.hmm_name, reconciled)
                searched.append(hit.hmm_name)

    # turn new_data into a new row (dataframe) to be concatenated
    new_row = pd.DataFrame(new_data, index=[genome_name])

    # concatenate new df with input df
    df = pd.concat([df, new_row])

    # Fill missing values with 0
    df.fillna(0, inplace=True)
    return df

# test populateDataFrame:
# df = pd.DataFrame()
# df = populateDataFrame(file_path, df)
# print(df)


"""
# Let's make a matrix (list of lists)
"""


def populateMatrix(file_path, matrix):
    """
    Parses the input faa file of Hmmer result, turns it into a "reconciled"
    dictionary of {gene: [HmmerHit]} using wholeGenomeOverlapRecon.

    Creates one "row" list of "domain frequencies" using the function domainSearch

    Returns a new matrix (list of lists) where the new "row" list is appended
    to the input matrix.
    """

    # create a dict {gene: [HmmerHit]} from file_path
    gene_hits = parseHmm(file_path)

    # extract genome name and add to matrix[1]
    pattern = '[A-Z]{2}_[A-Z]{3}_\d+\.\d+'
    matrix[0][re.search(pattern, file_path).group()] = len(matrix[0])

    # create a "reconciled" dict {gene: [HmmerHit]}
    reconciled = wholeGenomeOverlapRecon(gene_hits, 1e-5, 0.6)

    # Traverse reconciled hits; accumulate domainSearch results
    new_values = 20000 * [0]
    searched = []
    for gene, hits in reconciled.items():
        for hit in hits:
            # new domain
            if hit.hmm_name not in matrix[1]:
                matrix[1][hit.hmm_name] = len(matrix[1])
                new_values[len(matrix[1]) - 1] += domainSearch(hit.hmm_name, reconciled)
                searched.append(hit.hmm_name)
            # old domain
            elif hit.hmm_name not in searched:
                new_values[matrix[1][hit.hmm_name]] += domainSearch(hit.hmm_name, reconciled)
                searched.append(hit.hmm_name)

    # append new_values to the input matrix (list of lists)
    matrix.append(new_values)

    return matrix


# test populateMatrix for one file

# matrix = [{}, {}]
# matrix = populateMatrix(file_path, matrix)
# print(matrix)
# print(len(matrix[0]))


"""

========================================================
DS score
========================================================

"""
# make a dictionary of {phylum: [genomes]} from the ar53_taxonomy.tsv file
taxonomy_path = 'gtdb/ar53_taxonomy.tsv'
phylum_genomes = {}

with open(taxonomy_path, 'r') as file:
    output = file.readlines()

    for line in output:
        # Split the line into parts using tab ('\t') as the delimiter
        parts = line.split('\t')
        if len(parts) >= 2:
            # Extract the desired parts from the split line
            genome = parts[0]
            phylum = parts[1].split(';')[
                1]  # Assuming 'p__Phylum_name' is always the second part

            # Check if phylum is already in the dictionary
            if phylum in phylum_genomes:
                # Append the genome to the existing list
                phylum_genomes[phylum].append(genome)
            else:
                # Create a new list with the genome as the first element
                phylum_genomes[phylum] = [genome]


# how many phyla?
# phyla = list(phylum_genomes.keys())
# print(phyla)
# print(len(phyla))  # 20
#
# phylum_sizes = []
# for phylum, genomes in phylum_genomes.items():
#     phylum_sizes.append(len(genomes))
# print(phylum_sizes)  # [1429, 1220, 458, 1192, 62, 850, "11", 131, 137, 46, 229, 124, 77, "10", 20, 19, 16, 16, 6, 9]

# how many genomes?
# print(sum(phylum_sizes))  # 6062

# 3413 < 6062. What are we missing?
# matrix = pd.read_csv('domain_genome_matrix_archaea.csv', index_col=0)
# my_phylum_sizes = []
# for phylum, genomes in phylum_genomes.items():
#     num_genome = 0
#     for genome in genomes:
#         if genome in matrix.index:
#             num_genome += 1
#     my_phylum_sizes.append(num_genome)
# print(my_phylum_sizes)  # [768, 696, 172, 543, 32, 577, "1", 112, 92, 24, 176, 82, 58, "5", 16, 16, 14, 13, 6, 9]

"""
Create a big Dataframe 
"""
# df = pd.DataFrame()
# data_path = '/Users/tseamuscorlett/Desktop/LongoLab/Fold_Gated/protein_properties/gtdb/archaea'
# # change data_path and taxonomy_path above;
# # also change the csv file name below if needed
#
# # Iterate through the 'data' directory and process each .faa file
# for file_name in os.listdir(data_path):
#     faa_path = os.path.join(data_path, file_name)
#     df = populateDataFrame(faa_path, df)
#
# # Save the DataFrame to a CSV file
# df.to_csv('domain_genome_matrix_archaea.csv', index=False)

# DS (distribution score) of a domain = fraction of phyla for which 50% of the genomes in
# that phylum have the given domain (once for Archaea, once for Bacteria)

"""
Create a big matrix
"""
data_path = '/Users/tseamuscorlett/Desktop/LongoLab/Fold_Gated/protein_properties/gtdb/archaea'
matrix = [{}, {}]

# Iterate through the 'data' directory and process each .faa file
for file_name in os.listdir(data_path):
    faa_path = os.path.join(data_path, file_name)
    matrix = populateMatrix(faa_path, matrix)

# save matrix as csv
csv_file = "matrix_archaea.csv"

# Open the CSV file in write mode
with open(csv_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    for row in matrix:
        writer.writerow(row)



def calculateDS(domain_name, phyla, df, phylum_sizes=None):
    """
    Returns a DS (distribution score) of a domain_name (str) within the
    phyla {phylum: [genomes]}, where DS = fraction of phyla for which 50% of
    the genomes in that phylum have the given domain. Requires a dataframe
    as input where row_names = genome names, col_names = domain names, and values
    = frequency of domain in that genome

    Use custom list of phylum_sizes in case df does not contain all genomes
    """

    phylum_hit = 0
    counter = 0
    for phylum, genomes in phyla.items():
        genome_hit = 0
        for genome in genomes:
            # use Dataframe: Extract a value based on [column, row]
            if genome in df.index:
                if df.loc[genome, domain_name] > 0:
                    genome_hit += 1
        if phylum_sizes is None:  # no custom phylum_sizes
            if (genome_hit / len(phylum_genomes[phylum])) >= 0.5:
                phylum_hit += 1
        else:  # with custom phylum_sizes
            if (genome_hit / phylum_sizes[counter]) >= 0.5:
                phylum_hit += 1
            counter += 1
    return phylum_hit / len(phyla)


"""
test calculateDS for one domain:
"""
# matrix = pd.read_csv('domain_genome_matrix_archaea.csv', index_col=0)
# print(calculateDS('MCM_AAA', phylum_genomes, matrix))
# print(calculateDS('Sigma54_activat', phylum_genomes, matrix))


def calculateDSExtra(domain_name, phyla, df):
    """
    Returns a list [DS, [lost phyla], [lost genomes]]
    """

    lost_phyla = []
    lost_genomes = []
    phylum_hit = 0
    for phylum, genomes in phyla.items():
        genome_hit = 0
        for genome in genomes:
            # use Dataframe: Extract a value based on [column, row]
            if genome in df.index:
                if df.loc[genome, domain_name] > 0:
                    genome_hit += 1
                else:
                    lost_genomes.append(genome)
        if (genome_hit / len(phylum_genomes[phylum])) >= 0.5:
            phylum_hit += 1
        else:
            lost_phyla.append(phylum)
    return [phylum_hit / len(phyla), lost_phyla, lost_genomes]

"""
Figuring out why DS of Sigma54_activat = 0.9
by using calculateDSExtra
"""
# matrix = pd.read_csv('domain_genome_matrix_archaea.csv', index_col=0)

# print(calculateDSExtra('Sigma54_activat', phylum_genomes, matrix))
# [0.9, ['p__Huberarchaeota' (n=11), 'p__Undinarchaeota' (n=10)], ['GB_GCA_016784245.1', 'GB_GCA_016191585.1', 'GB_GCA_011047985.1', 'GB_GCA_002687735.1', 'GB_GCA_016186855.1']]
# where [p__Nanoarchaeota, p__Aenigmatarchaeota, no hit, p__Iainarchaeota, p__Iainarchaeota]

# print(calculateDSExtra('MCM_AAA', phylum_genomes, matrix))
# [0.9, ['p__Huberarchaeota' (n=11), 'p__Undinarchaeota' (n=10)], ['GB_GCA_016784245.1', ... LONG list]

# using the custom fourth argument in calculateDS to specify my_phylum_sizes
# my_phylum_sizes = [768, 696, 172, 543, 32, 577, 1, 112, 92, 24, 176, 82, 58, 5, 16, 16, 14, 13, 6, 9]
# print(calculateDS('Sigma54_activat', phylum_genomes, matrix, my_phylum_sizes))
# 1.0

# my_phylum_sizes = [768, 696, 172, 543, 32, 577, 1, 112, 92, 24, 176, 82, 58, 5, 16, 16, 14, 13, 6, 9]
# print(calculateDS('MoxR', phylum_genomes, matrix))  # 0.25
# print(calculateDS('MoxR', phylum_genomes, matrix, my_phylum_sizes))  # 0.55

# print(calculateDS('PPK2', phylum_genomes, matrix))  # 0.0
# print(calculateDS('PPK2', phylum_genomes, matrix, my_phylum_sizes))  # 0.0


"""
calculateDS for all domains:
"""
# DS_archaea_dict = {}
# matrix = pd.read_csv('domain_genome_matrix_archaea.csv', index_col=0)
# my_phylum_sizes = [768, 696, 172, 543, 32, 577, 1, 112, 92, 24, 176, 82, 58, 5, 16, 16, 14, 13, 6, 9]

# for domain_name in matrix.columns:
#     DS_archaea_dict[domain_name] = [calculateDS(domain_name, phylum_genomes, matrix, my_phylum_sizes)]
#
# # CSV file export
# DS_archaea = pd.DataFrame(DS_archaea_dict)
# DS_archaea.to_csv('DS_archaea.csv')
#
# DS_archaea = pd.read_csv('DS_archaea.csv', index_col=0)
# print(DS_archaea['MCM_AAA'])
# print(DS_archaea['Sigma54_activat'])
# print(DS_archaea['MoxR'])
# print(DS_archaea['fake domain name'])  # will raise KeyError


"""
How many domains in 3,413 genomes?
"""


def countDomains(file_path, domains):
    """

    """
    # create a dict {gene: [HmmerHit]} from file_path
    gene_hits = parseHmm(file_path)

    # iterate through genome files, populate "domains" list
    for gene, hits in gene_hits.items():
        for hit in hits:
            if hit.hmm_name not in domains.keys():
                domains[hit.hmm_name] = 1
            else:
                domains[hit.hmm_name] = domains[hit.hmm_name] + 1
    return domains

# data_path = '/Users/tseamuscorlett/Desktop/LongoLab/Fold_Gated/protein_properties/gtdb/archaea'
# domains = {}
# for file_name in os.listdir(data_path):
#     faa_path = os.path.join(data_path, file_name)
#     domains = countDomains(faa_path, domains)
# print(domains)
# print(len(domains))  # 10948


"""
Scaling law for each domain:
"""
# matrix = pd.read_csv('domain_genome_matrix_archaea.csv', index_col=0)
#
# # X-axis: make a list of genome sizes
# genome_sizes = []
# for index, row in matrix.iterrows():
#     # Calculate the sum of values in the current row
#     row_sum = row.sum()
#     genome_sizes.append(row_sum)
#
# # Y-axis: for one domain 'MCM_AAA'
# dom_freq_MCM_AAA = matrix['MCM_AAA'].tolist()
# # print(dom_freq_MCM_AAA)
#
# # Y-axis: for one domain 'Sigma54_activat'
# dom_freq_Sigma54_activat = matrix['Sigma54_activat'].tolist()
# # print(dom_freq_Sigma54_activat)
#
# # plot
# plt.scatter(genome_sizes, dom_freq_Sigma54_activat, marker='o', color='b', label='Data Points')
# plt.xlabel('genome size')
# plt.ylabel('domain count')
# plt.show()
# plt.savefig("images/Sigma54_activat.png")
