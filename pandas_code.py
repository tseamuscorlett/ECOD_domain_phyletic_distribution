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

"""
This file contains code using pandas dataframe
Removed from main after switching to using list of lists exclusively
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
original pupulateMatrix() using domainSearch()
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


def calculateDS(domain_name, phyla, df, phylum_sizes=None):
    """
    Returns a DS (distribution score) of an F/X-group (str) within the
    phyla {phylum: [genomes]}, where DS = fraction of phyla for which 50% of
    the genomes in that phylum have the given F/X-group. Requires a dataframe
    as input where row_names = genome names, col_names = domain names, and values
    = frequency of domain in that genome

    Use custom list of phylum_sizes in case df does not contain all genomes in phyla
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


def countPhylumSizes(matrix):
    my_phylum_sizes = []
    for phylum, genomes in phylum_genomes.items():
        num_genome = 0
        for genome in genomes:
            if genome in matrix.index:
                num_genome += 1
        my_phylum_sizes.append(num_genome)
    print(my_phylum_sizes)

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
# # CSV file export (pandas)
# DS_archaea = pd.DataFrame(DS_archaea_dict)
# DS_archaea.to_csv('DS_archaea.csv')
#
# DS_archaea = pd.read_csv('DS_archaea.csv', index_col=0)
# print(DS_archaea['MCM_AAA'])
# print(DS_archaea['Sigma54_activat'])
# print(DS_archaea['MoxR'])
# print(DS_archaea['fake domain name'])  # will raise KeyError

"""
map FtoX and save a new df as csv
"""

# archaea
# # use mapFtoX to create "x-group_genome_matrix_archaea.csv"
# df = pd.read_csv('domain_genome_matrix_archaea.csv', index_col=0)
#
# df = mapFtoX(df, f_to_x, same)
#
# # Save the DataFrame to a CSV file
# df.to_csv('x-group_genome_matrix_archaea.csv')


"""
DS for X-groups!
"""
# archaea
# my_phylum_sizes = [768, 696, 172, 543, 32, 577, 1, 112, 92, 24, 176, 82, 58, 5, 16, 16, 14, 13, 6, 9]
#
# x_group_DS_archaea_dict = {}
# x_matrix = pd.read_csv('x-group_genome_matrix_archaea.csv', index_col=0)
# for x_group in x_matrix.columns:
#     x_group_DS_archaea_dict[x_group] = [calculateDS(x_group, phylum_genomes, x_matrix, my_phylum_sizes)]
#
# # CSV file export (pandas)
# x_DS_archaea = pd.DataFrame(x_group_DS_archaea_dict)
# x_DS_archaea.to_csv('DS_archaea_xgroup.csv')


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