# Import Packages First
import csv
import math
import os
import scipy

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
    Requires first input matrix = [{}, {}]

    Parses the input faa file of Hmmer result, turns it into a "reconciled"
    dictionary of {gene: [HmmerHit]} using wholeGenomeOverlapRecon.

    Creates one "row" list of "domain frequencies" using the function domainSearch

    Returns a new matrix (list of lists) where the new "row" list is appended
    to the input matrix.
    """

    # create a dict {gene: [HmmerHit]} from file_path
    gene_hits = parseHmm(file_path)

    # extract genome name and add to matrix[0]
    pattern = '[A-Z]{2}_[A-Z]{3}_\d+\.\d+'
    matrix[0][re.search(pattern, file_path).group()] = len(matrix[0])

    # create a "reconciled" dict {gene: [HmmerHit]}
    reconciled = wholeGenomeOverlapRecon(gene_hits, 1e-5, 0.6)

    # Traverse reconciled hits; accumulate domainSearch results
    new_values = 10000 * [0]
    searched = set()
    for gene, hits in reconciled.items():
        for hit in hits:
            # new domain
            if hit.hmm_name not in matrix[1]:
                new_index = len(matrix[1])
                matrix[1][hit.hmm_name] = new_index
                new_values[new_index - 1] += domainSearch(hit.hmm_name, reconciled)
                searched.add(hit.hmm_name)
            # old domain
            elif not {hit.hmm_name}.isdisjoint(searched):
                new_values[matrix[1][hit.hmm_name]] += domainSearch(hit.hmm_name, reconciled)
                searched.add(hit.hmm_name)

    # append new_values to the input matrix (list of lists)
    matrix.append(new_values)

    return matrix


# test populateMatrix for one file

# matrix = [{}, {}]
# matrix = populateMatrix(file_path, matrix)
# print(matrix)
# print(len(matrix[0]))


"""
Let's make a matrix WITHOUT using domainSearch.
This improved speed for archaea by 2x (30min –> 15min) 
"""


def populateMatrixFast(file_path, matrix):
    """
    Requires first input matrix = [{}, {}]

    Parses the input faa file of Hmmer result, turns it into a "reconciled"
    dictionary of {gene: [HmmerHit]} using wholeGenomeOverlapRecon.

    Creates one "row" list of "domain frequencies" by traversing the reconciled hits.

    Returns a new matrix (list of lists) where the new "row" list is appended
    to the input matrix.
    """

    # create a dict {gene: [HmmerHit]} from file_path
    gene_hits = parseHmm(file_path)

    # extract genome name and add to matrix[0] 'dict of genomes'
    pattern = '[A-Z]{2}_[A-Z]{3}_\d+\.\d+'
    matrix[0][re.search(pattern, file_path).group()] = len(matrix[0])

    # create a "reconciled" dict {gene: [HmmerHit]}
    reconciled = wholeGenomeOverlapRecon(gene_hits, 1e-5, 0.6)

    # Traverse reconciled hits to update matrix[1] 'dict of domains'
    for gene, hits in reconciled.items():
        for hit in hits:
            if hit.hmm_name not in matrix[1]:
                matrix[1][hit.hmm_name] = len(matrix[1])

    # Traverse reconciled hits to accumulate 'hits per domain'
    new_values = 10000 * [0]
    for gene, hits in reconciled.items():
        for hit in hits:
            new_values[matrix[1][hit.hmm_name]] += 1

    # append new_values to the input matrix (list of lists)
    matrix.append(new_values)

    return matrix


# test populateMatrixFast for one file

# matrix = [{}, {}]
# matrix = populateMatrixFast(file_path, matrix)
# print(matrix)
# print(len(matrix[1]))


def parseMatrix(file_path):
    """
    Parse the output of populateMatrixFast
    Assumes the 1st row is the index row, the 2nd row is the column row
    """

    df = pd.read_csv(file_path, header=None)

    index_row = df.iloc[0]
    column_row = df.iloc[1]
    df = df[2:]  # df without the first two lines

    column_row.dropna(axis=0, inplace=True)  # drop NaN columns;
    index_row.dropna(axis=0, inplace=True)  # drop NaN rows;

    # Find the difference in df vs. column_row
    remove = len(df.columns) - len(column_row)
    # Cut the last n columns from df
    df = df.iloc[:, :-remove]

    df.index = index_row
    df.columns = column_row

    return df


def calculateDS(domain_name, phyla, df, phylum_sizes=None):
    """
    Returns a DS (distribution score) of an F/X-group (str) within the
    phyla {phylum: [genomes]}, where DS = fraction of phyla for which 50% of
    the genomes in that phylum have the given F/X-group. Requires a dataframe
    as input where row_names = genome names, col_names = domain names, and values
    = frequency of domain in that genome

    Use custom list of phylum_sizes in case df does not contain all genomes
    """

    phylum_hit = 0
    counter = 0
    for phylum, genomes in phyla.items():
        genome_hit = 0

        # search all genomes within phylum
        for genome in genomes:
            if genome in df.index:
                if df.loc[genome, domain_name] > 0:
                    genome_hit += 1

        # Does >50% of genomes in this phylum contain the domain?
        if phylum_sizes is None:  # no custom phylum_sizes
            if (genome_hit / len(phylum_genomes[phylum])) >= 0.5:
                phylum_hit += 1
        else:  # with custom phylum_sizes
            if (genome_hit / phylum_sizes[counter]) >= 0.5:
                phylum_hit += 1
            counter += 1
    return phylum_hit / len(phyla)


# test calculateDS for one domain:
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


def calculateDSMatrix(file_path, phyla, phylum_sizes=None):
    """
    Returns a dictionary of F/X-groups (str) : distribution_score (int)  within
    the phyla {phylum: [genomes]}, where DS = fraction of phyla for which 50% of
    the genomes in that phylum have the given domain.

    Requires a file_path of .csv file from which a MATRIX (list of lists) is made
    where 1st row is the genomes, 2nd row is the domain_names, and the rest are the
    domain frequencies.

    Use custom list of phylum_sizes in case MATRIX does not contain all genomes
    """

    with open(file_path, 'r') as file:
        output = file.readlines()

    # extract genomes_row & domains_row
    # remove ''s and the last '\n'
    genomes_row = list(filter(lambda x: (x != ''), output[0].split(',')))[:-1]
    domains_row = list(filter(lambda x: (x != ''), output[1].split(',')))[:-1]

    # turn lists into dicts for faster lookup
    genomes_dict = {item: index for index, item in enumerate(genomes_row)}
    domains_dict = {item: index for index, item in enumerate(domains_row)}

    # remove them from output
    output = output[2:]

    domainDS_dict = {}
    for domain in domains_dict:
        dom_index = domains_dict[domain]
        phylum_hit = 0
        counter = 0
        for phylum, genomes in phyla.items():
            genome_hit = 0

            # search all genomes within phylum
            for genome in genomes:
                if genome in genomes_dict.keys():
                    gen_index = genomes_dict[genome]
                    if int(output[gen_index].split(',')[dom_index]) > 0:
                        genome_hit += 1

            # Does >50% of genomes in this phylum contain the domain?
            if phylum_sizes is None:  # no custom phylum_sizes
                if (genome_hit / len(phylum_genomes[phylum])) >= 0.5:
                    phylum_hit += 1
            else:  # with custom phylum_sizes
                if (genome_hit / phylum_sizes[counter]) >= 0.5:
                    phylum_hit += 1
                counter += 1

        domainDS_dict[domain] = [phylum_hit / len(phyla)]
    return domainDS_dict


def countDomains(file_path, domains):
    """
    Updates and returns the given dict 'domains' by adding all domains found in
    the HMMER output file 'file_path'
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


def mapFtoX(df, f_to_x, same):
    """
    Takes df (rows = genomes, columns = F-groups) and returns a
    new_df (rows = genomes, columns = X-groups) using
    f_to_x {'xgroup': ['fgroup']} and same {'fgroup,fgroup': ['xgroup', 'xgroup']}
    """
    checked = {}
    no_match = []
    ambiguous = []

    for column in df.columns:
        # get X-group (f_to_x = {'F-group': ['X-group']})
        try:
            xgroup = f_to_x[column][0]  # e.g. '2004'
        # if no EXACT match (DALR_1,tRNA-synt_1d), check within 'same'
        except KeyError:
            for key in same:
                if column in key:
                    xgroup = same[key][0]
                    break  # some domains appear more than once... just take the first one
            else:
                print(f'no match for {column}')
                no_match.append(column)
                continue  # ignore this column

        # check that it's not "ambiguous"
        if xgroup != 'a':
            # new X-group
            if xgroup not in checked:
                checked[xgroup] = df[column].tolist()

            # X-group already exists => add the domain counts
            else:
                checked[xgroup] = [x + y for x, y in
                             zip(checked[xgroup], df[column])]
        if xgroup == 'a':
            ambiguous.append(column)
            print(f'{column} is ambiguous')

    new_df = pd.DataFrame(checked, index=[df.index])
    print(f'total # no match: {len(no_match)}')
    print(f'total # ambiguous F-group: {len(ambiguous)}')

    return new_df


def mapFtoXmatrix(file_path, f_to_x, same):
    """

    """

    with open(file_path, 'r') as file:
        output = file.readlines()

    # extract genomes_row & domains_row
    # remove ''s and the last '\n'
    genomes_row = list(filter(lambda x: (x != ''), output[0].split(',')))[:-1]
    domains_row = list(filter(lambda x: (x != ''), output[1].split(',')))[:-1]

    # turn lists into dicts for faster lookup
    genomes_dict = {item: index for index, item in enumerate(genomes_row)}
    domains_dict = {item: index for index, item in enumerate(domains_row)}

    # remove them from output
    output = output[2:]

    checked = {}
    no_match = []
    ambiguous = []

    for domain in domains_dict:
        # get X-group (f_to_x = {'F-group': ['X-group']})
        try:
            xgroup = f_to_x[domain][0]  # e.g. '2004'
        # if no EXACT match (DALR_1,tRNA-synt_1d), check within 'same'
        except KeyError:
            for key in same:
                if domain in key:
                    xgroup = same[key][0]
                    break  # some domains appear more than once... just take the first one
            else:
                print(f'no match for {domain}')
                no_match.append(domain)
                continue  # ignore this column

        # check that it's not "ambiguous"
        if xgroup != 'a':
            dom_index = domains_dict[domain][0]

            # new X-group
            if xgroup not in checked:
                checked[xgroup] = [x[dom_index] for x in output]

            # X-group already exists => add the domain counts
            else:
                checked[xgroup] = [x + y for x, y in
                                   zip(checked[xgroup], [x[dom_index] for x in output])]
        if xgroup == 'a':
            ambiguous.append(domain)
            print(f'{domain} is ambiguous')

    # save 'checked' dict as csv
    csv_file_path = "x-group_genome_matrix_archaea2.csv"

    # Write the dictionary to a CSV file
    with open(csv_file_path, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)

        # keep the genome names
        csv_writer.writerow(genomes_dict.keys())

        # Write the header (x-groups)
        csv_writer.writerow(checked.keys())

        # Write the data (f-groups)
        csv_writer.writerows(zip(*checked.values()))

    print(f'total # no match: {len(no_match)}')
    print(f'total # ambiguous F-group: {len(ambiguous)}')

    return

"""
Figuring out why DS of Sigma54_activat = 0.9
using calculateDSExtra
"""
# matrix = pd.read_csv('domain_genome_matrix_archaea.csv', index_col=0)

# print(calculateDSExtra('Sigma54_activat', phylum_genomes, matrix))
# [0.9, ['p__Huberarchaeota' (n=11), 'p__Undinarchaeota' (n=10)], ['GB_GCA_016784245.1', 'GB_GCA_016191585.1', 'GB_GCA_011047985.1', 'GB_GCA_002687735.1', 'GB_GCA_016186855.1']]

# print(calculateDSExtra('MCM_AAA', phylum_genomes, matrix))
# [0.9, ['p__Huberarchaeota' (n=11), 'p__Undinarchaeota' (n=10)], ['GB_GCA_016784245.1', ... LONG list]

# using the custom fourth argument in calculateDS to specify my_phylum_sizes
# my_phylum_sizes = [768, 696, 172, 543, 32, 577, 1, 112, 92, 24, 176, 82, 58, 5, 16, 16, 14, 13, 6, 9]

# print(calculateDS('Sigma54_activat', phylum_genomes, matrix, my_phylum_sizes))  # 1.0

# print(calculateDS('MoxR', phylum_genomes, matrix))  # 0.25
# print(calculateDS('MoxR', phylum_genomes, matrix, my_phylum_sizes))  # 0.55

# print(calculateDS('PPK2', phylum_genomes, matrix))  # 0.0
# print(calculateDS('PPK2', phylum_genomes, matrix, my_phylum_sizes))  # 0.0


"""
========================================================
DS score
========================================================
"""
# make a dictionary of {phylum: [genomes]} from taxonomy.tsv file
taxonomy_path = 'gtdb/ar53_taxonomy.tsv'
# taxonomy_path = 'gtdb/bac120_taxonomy.tsv'
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


# # how many phyla?
# phyla = list(phylum_genomes.keys())
# # print(phyla)
# print(len(phyla))  # 20 archaea, 169 bacteria

# phylum_sizes = []
# for phylum, genomes in phylum_genomes.items():
#     phylum_sizes.append(len(genomes))
# print(phylum_sizes)
# archaea: [1429, 1220, 458, 1192, 62, 850, "11", 131, 137, 46, 229, 124, 77, "10", 20, 19, 16, 16, 6, 9]
# bacteria: [141114, 61795, 28532, 21744, 6845, 2331, 20893, 2848, 612, 2818, 1230, 1472, 412, 295, 1503, 2214, 4645, 31, 101, 492, 208, 496, 499, 627, 42, 276, 456, 228, 313, 109, 120, 122, 1602, 147, 404, 244, 116, 22, 10, 325, 63, 21, 414, 79, 223, 171, 54, 54, 21, 119, 62, 102, 105, 92, 53, 107, 54, 161, 7, 12, 82, 24, 92, 32, 11, 77, 50, 52, 16, 10, 12, 14, 7, 55, 50, 10, 15, 25, 20, 14, 8, 8, 29, 25, 38, 7, 27, 27, 22, 11, 23, 15, 5, 25, 27, 8, 4, 17, 2, 6, 3, 2, 4, 7, 4, 5, 4, 2, 2, 10, 6, 3, 7, 3, 4, 4, 10, 3, 8, 2, 5, 7, 2, 2, 1, 6, 3, 8, 1, 2, 2, 2, 1, 11, 4, 4, 1, 3, 1, 3, 1, 2, 1, 2, 2, 3, 2, 3, 1, 3, 1, 2, 2, 1, 2, 1, 1, 1, 2, 5, 1, 1, 2, 2, 1, 1, 1, 1, 1]

# how many genomes?
# print(sum(phylum_sizes))  # 6062

# 3412 < 6062. What are we missing?
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
Create a big matrix (FAST)
"""
# data_path = '/Users/tseamuscorlett/Desktop/LongoLab/Fold_Gated/protein_properties/gtdb/archaea'
# matrix = [{}, {}]
#
# # Iterate through the 'data' directory and process each .faa file
# for file_name in os.listdir(data_path):
#     faa_path = os.path.join(data_path, file_name)
#     matrix = populateMatrixFast(faa_path, matrix)
#
# # save matrix as csv
# csv_file = "matrix_archaea_fast.csv"
#
# # Open the CSV file in write mode
# with open(csv_file, mode='w', newline='') as file:
#     writer = csv.writer(file)
#     for row in matrix:
#         writer.writerow(row)

"""
calculateDS for all domains:
"""
# archaea

# DS_archaea_dict = {}
# matrix = pd.read_csv('domain_genome_matrix_archaea.csv', index_col=0)
# my_phylum_sizes = [768, 696, 172, 543, 32, 577, 1, 112, 92, 24, 176, 82, 58, 5, 16, 16, 14, 13, 6, 9]
#
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


# # archaea with calculateDSMatrix (DON'T use pandas df)
# my_phylum_sizes = [768, 696, 172, 543, 32, 577, 1, 112, 92, 24, 176, 82, 58, 5, 16, 16, 14, 13, 6, 9]
# DS_archaea_dict = calculateDSMatrix('matrix_archaea_fast.csv', phylum_genomes, my_phylum_sizes)
#
# # Specify the CSV file path
# csv_file_path = 'output.csv'
#
# # Write the dictionary to a CSV file
# with open(csv_file_path, 'w', newline='') as csv_file:
#     csv_writer = csv.writer(csv_file)
#
#     # Write the header
#     csv_writer.writerow(DS_archaea_dict.keys())
#
#     # Write the data
#     csv_writer.writerows(zip(*DS_archaea_dict.values()))


# bacteria
# my_phylum_sizes = [get this from matrix]
# DS_bacteria_dict = {}

# not enough RAM!! Running this forces Python to halt the run:
# parsed_bac_matrix = parseMatrix('matrix_bacteria_fast.csv')
# parsed_bac_matrix.to_csv('matrix_bacteria_fast_parsed.csv')




"""
How many domains in archaea 3,412 genomes?
"""

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

"""
F –> X-groups
"""

# Create a dict f_to_x = {domain.f_name: ['X-group']} from ecod.develop279.domains.txt
# Make sure to record "ambiguous" case of 1 F-group mapping to 2+ X-groups

# with open('data/ecod.develop279.domains.txt', 'r') as file:
#     output = file.readlines()
# f_to_x = {}
# for line in output:
#     if line[0] == '#':  # skip the first lines
#         continue
#     domain = EcodDomain(line)
#     # new domain (F-group) name
#     if domain.f_name not in f_to_x:
#         f_to_x[domain.f_name] = [domain.f_id.split('.', 1)[0]]  # 2004.1.1.140 => 2004
#     # old domain (F-group) name AND the X-group is different
#     elif f_to_x[domain.f_name][0] != domain.f_id.split('.', 1)[0]:
#         f_to_x[domain.f_name] = 'ambiguous'  # updates the value to 'ambiguous'

# print(f_to_x)
# print(len(f_to_x))  # 12634: ECOD v.279

# # Write the dictionary to a CSV file
# csv_file_path = 'f_to_x.csv'
# with open(csv_file_path, 'w', newline='') as csv_file:
#     csv_writer = csv.writer(csv_file)
#
#     # Write the header
#     csv_writer.writerow(f_to_x.keys())
#
#     # Write the data
#     csv_writer.writerows(zip(*f_to_x.values()))

"""
check cases where f_name contains a comma (e.g., DALR_1,tRNA-synt_1d)
"""

# with open('data/ecod.develop279.domains.txt', 'r') as file:
#     output = file.readlines()
# commas = set()
# for line in output:
#     if line[0] == '#':  # skip the first lines
#         continue
#     domain = EcodDomain(line)
#     if ',' in domain.f_name:
#         commas.add(domain.f_name)
#
# print(commas)

# save the set
# commas = {'SecA_PP_bind,SBP_bac_1_C', 'putA_2nd,Pro_dh-DNA_bdg', 'KOG2546,ATP-synt_ab_Xtn', 'GIDA_assoc_1st,FAD_binding_3_1st', 'HD_5,RNAse_Pc_like', '7tm_2,Rubredoxin', 'GF_recep_IV,Recep_L_domain', 'XRCC4_3rd_2,Tropomyosin_2', 'RNAse_Pc_like,GATA', 'NCD3G,ANF_receptor_2nd_1', 'Hormone_recep_1,SBP_bac_1_N_1', 'NNMT_PNMT_TEMT,RNAse_Pc_like', 'BPD_transp_1,malF', 'PRK11000,ABC_tran_1', 'PLN02303_like,Urease_beta', 'Pyr_redox_dim,Pyr_redox_2', 'PYC_OADA_2nd,PYC_OADA_1st,HMGL-like_2', 'Neur_chan_LBD,Neur_chan_memb', 'PHD,BAH_1', 'Reo_sigma1_C,Reo_sigma1_N', 'PLN00104,RNAse_Pc_like', 'TNFRSF11A_2nd,TNFRSF21', 'zf-CCCH_4,PHD3_KDM5A_like', 'lambda-1,zf-C2H2_6', 'PRK15442_1,Beta-lactamase2', 'BolA,RNAse_Pc_like', 'KOG0384_1,Helicase_C', 'Mst1_SARAH,V-set', 'Internalin_N,LRR_8_3', 'SBP_bac_3_1st,SBP_bac_3_2nd', 'Viral_DNA_bp_2nd,Viral_DNA_bp_3rd', 'RNAse_Pc_like,Vinculin_1st', 'VGCC_beta4Aa_N,SH3_1_2', 'PSI_1,Sema_1', 'LbR_YadA-like,KOG0837_1', 'Sad1_UNC,EUF07923', 'Pkinase_Tyr,UBA_AID_AMPKalpha', 'PRK00566_3rd,RNA_pol_Rpb1_5_3rd', 'Gal_mutarotas_2,DUF5110_1,Glyco_hydro_31_N', 'Ank_2,Arm', 'GMC_oxred_C,GMC_oxred_N', 'PAD,PAD_M', 'KOG4642,TPR_11_4', 'DUF2094,RNAse_Pc_like', 'PLAT,RNAse_Pc_like', 'OTCace_N,RNAse_Pc_like', 'Pectate_lyase_3,Glyco_hydro_120', 'DNA_gyraseB,HATPase_c_1', 'tRNA-synt_1c_C_C,tRNA-synt_1c', 'tRNA-synt_1_1st,leuS_2nd', 'Beta-lactamase_1st,Beta-lactamase_2nd_1,RNAse_Pc_like', 'PRK00566_2nd,RNA_pol_Rpb1_5_3rd', 'SAM_Arap1,2,3_like', 'ANAPC4_WD40_9,Sec16_C', 'PSI_2,EUF07736', 'HD,PolyA_pol_RNAbd', 'KOG0196,Ephrin_lbd', 'PRK11000,ABC_tran', 'Ribonuc_red_lgN,ATP-cone_1', 'PPV_E2_N_C,PPV_E2_N_N', 'HypA,zf-C2H2_19', 'LAGLIDADG_1,Spindle_Spc25', 'MDA5_ID,Helicase_C_1', 'Cytochrome_C554,Cytochrom_C552', 'Methyltransf_11_3,RNAse_Pc_like', 'Hpt,Ank_2', 'Ion_trans_3,KOG0837_1', 'UDPGP_C,UDPGP_N', 'PRK15442,Peptidase_S11_1st_1', 'Hydrolase_4,Peptidase_S9_N', 'Aldo_ket_red,RNAse_Pc_like', 'AMP-binding_2nd,AMP-binding_3rd,AMP-binding_C', 'Fe-ADH_2,Fe-ADH_N', 'TPR_15,RNAse_Pc_like', 'Corona_nucleoca_1st,RNAse_Pc_like', 'I-set_9,I-set_13', 'BACK,BTB', 'murQ_2nd,PRK12570', 'RAI1,RNAse_Pc_like', 'LAGLIDADG_3,RNAse_Pc_like', 'Transferrin_2nd,Transferrin_1st', 'Gcd10p_1st,Gcd10p_2nd', 'RNAse_Pc_like,MarR', 'MR_MLE_N,RNAse_Pc_like', 'KOG2181,LIM_C_3', 'KOG0948,rRNA_proc-arch_3rd', 'GLF,Amino_oxidase_1st', 'ADH_zinc_N,ADH_N_2_1', 'Thrombin_light_1,Trypsin', 'Anticodon_1,tRNA-synt_1g', 'DNA_topoisoIV_2nd_1,DNA_topoisoIV_3rd', 'zf-C5HC2,KOG1246_1st', 'CdAMP_rec,RNAse_Pc_like', 'ANATO,A2M_N_2', 'PROCN_C,PROCN_N', 'PNPase,RNase_PH_1', 'Bac_luciferase,RNAse_Pc_like', 'MOSC_N,Phage_lysozyme_2', 'SAS-6_N_1,EUF08200', 'FOLN_like_1,EUF08727', 'DUF3443_like,RNAse_Pc_like', 'C2,KOG1265', 'PSI_3,EUF07736', 'Glyco_hydro_1,RNAse_Pc_like', 'Peptidase_S6_C,Glyco_hydro_120', 'Ribosomal_L19e_N,EUF07990', 'WIYLD,RNAse_Pc_like', 'Inhibitor_I29,RNAse_Pc_like', 'Metallophos_2_1,EUF08528', 'RRM_1,RNAse_Pc_like', 'EUF08255,SelA_N', 'Glycohydro_20b2,Glyco_hydro_20', 'KOG2181_1,LIM_C_3', 'YadA_stalk_1,EUF07855', 'RNAse_Pc_like,Clathrin_H_link', 'Molybdop_Fe4S4,Molybdopterin_2nd', 'KOG2292,RNAse_Pc_like', 'PksD,EUF08630', 'Molybdopterin_2nd,Molybdop_Fe4S4_1', 'DUF592,SIR2_1st', 'ATP-synt_ab,ATP-synt_ab_N', 'DAP3_C,DAP3_N', 'DUF1861,RNAse_Pc_like', 'EUF07844,Rieske_6', 'Peptidase_M18_1st,Peptidase_M20_1st', 'NOB1_Zn_bind,Ribosomal_L33', 'Glyco_hydro_66_1st,RNAse_Pc_like', 'Acetyltransf_10_4,SCP2_2_N', 'KOG2680_like,TIP49_3rd_1', 'Ribonuc_red_lgC,Ribonuc_red_lgN_1', 'KOG0947_3rd,KOG0948', 'DehI,RNAse_Pc_like', 'NUDIX_1,nudC_2nd', 'Adenosine_kin,RNAse_Pc_like', 'SelB-wing_1,SelB-wing_2', 'SBP_bac_1_C,SBP_bac_1_N_2', 'zf-CCHC,PknG_rubred_like', 'Fusion_gly_2,KOG0837_1', 'gyrB_1st_1,gyrB_1st_2', 'HA2_N,DEAD_3', 'Tudor-knot,RNAse_Pc_like', 'PNPase_1,RNase_PH_1', 'Gal_mutarotas_2,Glyco_hydro_31_N', 'E1-E2_ATPase_N,ATPase-IIC_X-K', 'EUF07859,Reo_sigma1_N_1', 'XRCC4_3rd_1,KOG0161_2nd', 'SIR2_2nd_1,SIR2_1st', 'EBP50_C,PDZ', 'Methyltransf_15,Met_10', 'KOG0652_3rd,KOG0652_2nd', 'RNAse_Pc_like,Enoyl_reductase', 'MGS,EUF08629', 'Pkinase_Tyr,RGS', 'PDZ_2,Trypsin_2', 'H2O2_YaaD,MRP-S27', 'Amidohydro_1_N,Amidohydro_1_C_2', 'KOG0275_2nd,zf_topless', 'Amino_oxidase_1st,Pyr_redox_2_2', 'F5_F8_type_C,RNAse_Pc_like', 'PG_binding_1,EUF07947', 'PRK09243,NAPRTase', 'RNAse_Pc_like,DisA_N', 'Amidohydro_3,D-HYD_N', 'KOG1246_1st,JmjC', 'DALR_1,tRNA-synt_1d', 'EUF08400,vigilin_like_KH_like', 'LRR_8_4,F-box-like_3', 'RNAse_Pc_like,OpuAC_1', 'Hemagglutinin_2nd,Hemagglutinin_1st,EUF08590', 'Sigma54_activat,CDC48', 'RNAse_Pc_like,Dimerisation2_like', 'Dimerisation2,RNAse_Pc_like', 'SBP_bac_1_N_1,RNAse_Pc_like', 'Adap_comp_sub_2nd,Adap_comp_sub_1st', 'PYC_OADA_1st,PYC_OADA_2nd_1', 'YciM,PHD_1', 'Fe_hyd_lg_C,Fe_hyd_SSU', 'KOG1155,ANAPC8', 'TPR_19,ANAPC8_1', 'Sigma70_r3,Sigma70_r2_2', 'TPR_19,ANAPC8', 'Asp_C,Asp_N', 'Methyltransf_11_4,Met_10', 'UvrA,ABC_tran', 'GalP_UDP_transf,GalP_UDP_tr_C', 'PA14,RNAse_Pc_like', 'GP41_2nd_1,GP41_like', 'PBP2_iGluR_NMDA_Nr2_like_1st,Lig_chan_5', 'Granulin_1,EUF08737', 'LepB_N,EUF08529', 'Raptor_N,KOG1517,Arm', 'TPK_catalytic,RNAse_Pc_like', 'DALR_2,tRNA-synt_1g', 'AMP-binding_2nd,AMP-binding_3rd', 'Myosin_head,IQ_2', 'PBP_dimer_2,Transpeptidase_2nd', 'PRK05762_1st,DNA_pol_B_exo1_3rd', 'B,Pkinase_Tyr', 'Flu_NP_1st,Flu_NP_2nd', 'ATP-synt_ab,ATP-synt_ab_C', 'Hemagglutinin_2nd,Hemagglutinin_1st_1', 'Lipase_chap,PolyA_pol_RNAbd', 'Calpain_III,RNAse_Pc_like', 'PRK04342_2nd,PRK04342_1st', 'EAL,RNAse_Pc_like', 'DUF4091,Glyco_hydro_123_N', 'tRNA-synt_1f_1st,LysS', 'HSP70_2nd,HSP70_3rd', 'LppC_1st,RNAse_Pc_like', 'GP41_2nd_1,EUF08590', 'KOG1038,DUF4909', 'ITAM_like,SH3_1_1', 'I-set_13,I-set', 'TusA,RNAse_Pc_like', 'NUDIX_1,KOG0648', 'UxuA,RNAse_Pc_like', 'Rieske_cytochrome_bc1,Rieske_6', 'EUF08311,RNAse_Pc_like', 'GntR,Transcrip_reg_1st', 'Flavi_glycoprot,Flavi_glycop_C', 'EF-hand_7_8,RNAse_Pc_like', 'TIG_4,HLH', 'Ribosomal_S13_N,Ribosomal_S15_2', 'Flo11,RNAse_Pc_like', 'zf-MYND,SET', 'KOG0947_3rd,DSHCT', 'GFP,zf-C2H2_4_like', 'Trp_DMAT,RNAse_Pc_like', 'EF-hand_7_8,IQ_1', 'RhbC_1,IucA_IucC', 'DUF1730,EUF08225', 'KOG1038,RPOL_N', 'RNA_pol_Rpb2_6_2nd,Sigma70_r4', 'DNA_pol3_beta,RNAse_Pc_like', 'Ribonuc_red_sm,RNAse_Pc_like', 'RNAse_Pc_like,HRM_1', 'Sortilin_C_1st,Sortilin-Vps10', 'DUF2225,KOG1086', 'Mfa_like_2_C,Mfa_like_2_N', 'Synuclein,SBP_bac_1_N_1', 'B,KOG0837_1', 'Raptor_N,HEAT_2', 'CODH,Prismane_2nd', 'LacAB_rpiB,RNAse_Pc_like', 'PhoD_1,SapB_1', 'STAT_bind_C,STAT_bind_N', 'NAD_binding_2,RNAse_Pc_like', '7tm_1_1,Cytochrom_B562_1', 'Glyoxalase_2_C,Glyoxalase_14', 'GHMP_kinases_N,RNAse_Pc_like', 'Integrin_beta,KOG1226_1st', 'KOG1562,Spermine_synt_N_1', 'Rho_N,Rho_RNA_bind', 'Secretin_N,PRK15339', 'Pterin_bind,RNAse_Pc_like', 'PNPase_KH,RNase_PH_C_2', 'Chlam_vir_2nd,Chlam_vir_3rd', 'Flu_PB2_2nd,Flu_PB2_3rd', 'dUTPase,RNAse_Pc_like', 'Calici_coat_C_N,Calici_coat_C_C', 'ERM,FERM_M', 'Lig_chan,Lig_chan_1', 'SdrD_B,RNAse_Pc_like', 'DALR_1,tRNA-synt_1_1st', 'Terminase_5,EUF08177', 'RNAse_Pc_like,C2_1', 'Amidohydro_1_C_2,PRK09229', 'Thiol_cytolysin,RNAse_Pc_like', 'Molybdop_Fe4S4,arsenite_ox_L', 'Longin,Synaptobrevin', 'RNAse_Pc_like,AhpC-TSA_1', 'Spectrin,EUF08159', 'zf-UDP,Rep_fac-A_C_N', 'OsmC,RNAse_Pc_like', 'RsmF_methylt_CI,EUF07530', 'I-set_4,fn3', 'ANF_receptor_2nd,ANF_receptor_1st', 'DEAD_1,KOG0951_2nd', 'Glyco_hydro_20b,RNAse_Pc_like', 'HAMP,HisKA', 'Sec23_BS,Sec23_trunk', 'PHF5_1,PHF5_2', 'Glyco_hydro_56,hEGF_3', 'LpxD,RNAse_Pc_like', 'Peptidase_S13_1st,Peptidase_S13_2nd', 'Kringle,WSC', 'gyrB_2nd,gyrB_1st_1,gyrB_1st_2', 'SBP_bac_1_N_1,PRK15442_1', 'Sacchrp_dh_C,Sacchrp_dh_NADP', 'Chromo,RNAse_Pc_like', 'Pox_polyA_pol_C,Pox_polyA_pol', 'Gla_like,Gla', 'Reo_sigma1_C,EUF07859,Reo_sigma1_N_1', 'PRORP_N,PRORP_C', 'RNAse_Pc_like,YcgR_like', 'BNR_2,7tm_1_1', 'Complex1_49kDa_1st,NiFeSe_Hases_1st', 'Crystallin,HSP20', 'PUF,Spermine_synth', 'ketoacyl-synt,Docking_1', 'AHSA1,RNAse_Pc_like', 'p450,RNAse_Pc_like', 'KOG1464,PCI_1', 'Fusion_gly_2,KOG0837', 'TNFR_c6_2,TNFRSF5', 'KOG0196,Laminin_EGF', 'Inhibitor_I69,Peptidase_C10', 'HA2_N,HA2_C', 'PP2C,RNAse_Pc_like', 'PYC_OADA_2nd,PYC_OADA_1st', 'Lyase_aromatic_C,Lyase_aromatic_N', 'OmpH,Vinculin_2nd', 'Fer4_7_1,NADH-G_4Fe-4S_3', 'P-II,PGM_PMM_IV', 'zf-C2H2_1,HSP33_C', 'MerR_1,TipAS', 'DUF5110,Glyco_hydro_31_C', 'Glycophorin_A,SSF', 'PQQ_like,Beta-lactamase2', 'ApbA,ApbA_C', 'ChAPs_C,ChAPs_N', 'I-set_15,I-set_11', 'PRK10954,Filamin_1', 'flgG,Flg_bb_rod', 'Raptor_N,Arm', 'UPF0146,RNAse_Pc_like', 'MukE_C,MukE_N', 'AAA_33,ADK_lid', 'Met_10,Methyltransf_11_2', 'hEGF_3,AMA-1_2nd', 'KOG0161_like,XRCC4_3rd_1', 'Anticodon_1,tRNA-synt_1_1st', 'PADR1,PLN03123', 'Pyridox_ox_2,Yop-YscD_ppl_2nd', 'TMAO_torS,EUF08144', 'ThiG,TIM', 'RNAse_Pc_like,Acyl_transf_1_1st', 'GDP_Man_Dehyd_1,RNAse_Pc_like', 'Ank_2,RNAse_Pc_like', 'ADH_zinc_N,ADH_N', '34_2nd,34_1st', 'B,Phage_lysozyme_2', 'LRR_8_4,KOG4308_like', 'PHD_BAZ1A_like_like,PSI', 'PRK09198,NAPRTase', 'Cytochrom_B561,YajC', 'TPP_enzyme_N,RNAse_Pc_like', 'Translin,RNAse_Pc_like', 'murQ_2nd,murQ_1st', 'Torsin,RecA', 'Glyoxalase_17,Glyoxalase_11', 'FKBP_C,RNAse_Pc_like', 'HTH_20_like,Lactamase_B_1', 'ACT_2,RNAse_Pc_like', 'GalA,RNAse_Pc_like', 'T6SS_VipA,VipB_N', 'Caudo_bapla_RBP_1st_1,Caudo_bapla_RBP_1st', 'Sigma54_activat,Vps4_C_C', 'Cytidylate_kin_2,Thymidylate_kin', 'PSI_2,TIG_3', 'Rieske,HcaE', 'EF-hand_7_8,IQ_like', 'Sigma54_activat,KOG0729', 'YfaS_3rd,YfaS_4th', 'RcbX,RNAse_Pc_like', 'KOG4008,Peripla_BP_4_C', 'DnaJ,RNAse_Pc_like', 'ThiF_2,APPBP1_RUB', 'Nucleoporin_C_C,Nucleoporin_C_N', 'LRR_8_4,F-box-like_4', 'CARD,RNAse_Pc_like', 'EF-hand_7_1,TFA2_2nd', 'RNAse_Pc_like,Hsm3_like', 'Methyltransf_11_2,Gcd10p_2nd', 'Beta-lactamase_1st,Beta-lactamase_2nd_1', '2-Hacid_dh,DUF3410', 'Integrin_b_cyt,IRS_1', 'Sec23_BS,zf-Sec23_Sec24', 'PksD,Acyl_transf_1_1st', 'RNAse_Pc_like,Lactamase_B_6_1', 'EUF08675,NAPRTase', 'Beta-prism_lec,RNAse_Pc_like', 'EUF07924,EUF07925', 'IPPT,IPT', 'RNA_pol_Rpb2_7,RNA_pol_Rpb1_5_3rd_1,RNA_pol_Rpb1_1', 'SET_5,zf-MYND', 'DAO_1st,FMO-like_2nd', 'Pentapeptide_4,RNAse_Pc_like', 'Ribosomal_S5_3,Ribosomal_S5_C_1', 'Aminotran_3_C_1,Aminotran_3_N', 'Thrombin_light,Trypsin', 'EUF08456,RNAse_Pc_like', 'LRR_8_3,RHH_6', 'CTP_transf_like_1,Lipoprotein_6', 'KOG4405,GDI_1st', 'Methyltransf_11_1,RNAse_Pc_like', 'TAXi_C,Asp_N', 'zf-C2H2_18,zf-H2C2_2_like', 'HeLo,RNAse_Pc_like', 'TPR_11_6,RNAse_Pc_like', 'Gal_mutarotas_2,Trefoil', 'RNAse_Pc_like,adh_short_C2_1', 'BLM10_mid,CLASP_N_2', 'DUF892,EUF08144', 'tRNA-synt_1g,Anticodon_Ia_Cys_like', 'PRK15442_1,Peptidase_S11_1st_1', 'PRK15442,Beta-lactamase2', 'KOG1226_1st,hEGF_like_1', 'KOG0163_1,Myosin_head_1', 'Dirigent,EUF08472', 'AAA_33,RNA_pol_Rpc34', 'fn3,EUF08472', 'PYC_OADA_2nd,HMGL-like_2', 'KR_1_SDR_x,EUF08195', 'Bunya_RdRp_3rd,EUF08583', 'Alpha_E2_glycop_2nd,Alpha_E2_glycop_3rd', 'DNA_RNApol_7kD,Shisa', 'AMP-binding_1st,AMP-binding_2nd', 'ZF_C2H2,Rep_fac-A_C_N', 'GFP,DUF4414', 'Plexin_cytopl_1st,KOG0837_1', 'zf-C2H2_13,zf-H2C2_2_like', 'Pico_P1A_C,Pico_P1A_N', 'ADK_lid,Cytidylate_kin_2', 'PCI,KOG2581', 'Se-cys_synth_N_like,EUF08255', 'SATase_N,Hexapep_2_1', 'DUF4982,Lectin_C', 'C1-set_3,C2-set_2_1', 'Peptidase_M27,RNAse_Pc_like', 'CDO_I,RNAse_Pc_like', 'BLM10_mid,BLM10_N', 'SH2,RNAse_Pc_like', 'Ank_2,FERM_f0', 'B,GP41_2nd', 'GGDEF,KOG0837_1'}

# # let's split them
# xpairs = []
# for comma in commas:
#     xpair = comma.split(',', 1)
#     # f_name contains 1 comma
#     if ',' not in xpair[1]:
#         xpairs.append(xpair)
#     # f_name contains 2+ commas
#     else:
#         last_two = xpair[1].split(',', 1)
#         new = [xpair[0], *last_two]  # unlist using *
#         xpairs.append(new)
# print(xpairs)

# # find their X-group integer code to see if they're the same or not
# xpairs_int = {}
# for pair in xpairs:
#     integers = []
#     for domain in pair:
#         # search within f_to_x dictionary
#         try:
#             integers.append(f_to_x[domain][0])
#         # domain name NEVER appears in isolation (e.g., KOG1038)
#         except:
#             for key in f_to_x:
#                 if domain in key:
#                     integers.append(f_to_x[key][0])
#                     print(key)
#                     break  # only need one
#     # remove anything with 'ambiguous'
#     if 'a' not in integers:
#         if len(integers) == 3:
#             xpairs_int[pair[0] + ',' + pair[1] + ',' + pair[2]] = integers
#         if len(integers) == 2:
#             xpairs_int[pair[0] + ',' + pair[1]] = integers
# print(xpairs_int)
# print(len(xpairs_int))  # 313 (436 including pairs with 'a')
#
# # create a dictionary of double/triple F_name entries with the SAME X-group
# same = {}
# for xpair, domains in xpairs_int.items():
#     if len(domains) == 2:
#         if domains[0] == domains[1]:
#             same[xpair] = domains
#     if len(domains) == 3:
#         if domains[0] == domains[1] == domains[2]:
#             same[xpair] = domains

# print(len(same))  # 87
# print(same)

# Write the dictionary to a CSV file
# csv_file_path = 'same.csv'
# with open(csv_file_path, 'w', newline='') as csv_file:
#     csv_writer = csv.writer(csv_file)
#
#     # Write the header
#     csv_writer.writerow(same.keys())
#
#     # Write the data
#     csv_writer.writerows(zip(*same.values()))


"""
Finally, let's map FtoX and save a new matrix as csv
"""

# # let's use mapFtoX to create "x-group_genome_matrix_archaea.csv"
# df = pd.read_csv('domain_genome_matrix_archaea.csv', index_col=0)
#
# df = mapFtoX(df, f_to_x, same)
#
# # Save the DataFrame to a CSV file
# df.to_csv('x-group_genome_matrix_archaea.csv')


"""
DS for X-groups!
"""
# archaea (make sure phylum_genomes comes from archaea_taxonomy.tsv)
# my_phylum_sizes = [768, 696, 172, 543, 32, 577, 1, 112, 92, 24, 176, 82, 58, 5, 16, 16, 14, 13, 6, 9]
#
# x_group_DS_archaea_dict = {}
# x_matrix = pd.read_csv('x-group_genome_matrix_archaea.csv', index_col=0)
# for x_group in x_matrix.columns:
#     x_group_DS_archaea_dict[x_group] = [calculateDS(x_group, phylum_genomes, x_matrix, my_phylum_sizes)]
#
# # CSV file export
# x_DS_archaea = pd.DataFrame(x_group_DS_archaea_dict)
# x_DS_archaea.to_csv('DS_archaea_xgroup.csv')



# bacteria (make sure phylum_genomes comes from bacteria_taxonomy.tsv)
# my_phylum_sizes = [get this from matrix]
# x_DS_bacteria_dict = calculateDSMatrix()
