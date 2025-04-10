import csv

def dict2csv(dict, csv_file_path):
    # Write the dictionary to a CSV file
    with open(csv_file_path, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)

        # Write the header
        csv_writer.writerow(dict.keys())

        # Write the data
        csv_writer.writerow(dict.values())


def csv2dict(csv_file_path):
    with open(csv_file_path, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)

        # Assuming the first row contains the headers
        headers = next(csv_reader)

        # Initialize an empty dictionary with each header as a key and an empty list as the value
        result_dict = {header: [] for header in headers}

        # Iterate through the rows and append values to the corresponding lists
        for row in csv_reader:
            for header, value in zip(headers, row):
                result_dict[header].append(value)

    return result_dict


def parseMatrix(file_path):
    """
    Parse the matrix (list of lists) of domain vs. genome HMMER-hit data
    Assumes the 1st row is the index row, the 2nd row is the column row

    Returns domains_dict, genomes_dict, output (list of lists)
    """

    with open(file_path, 'r') as file:
        output = file.readlines()

    # extract genomes_row & domains_row
    # remove ''s and the last '\n'
    genomes_row = list(filter(lambda x: (x != ''), output[0].strip().split(',')))
    domains_row = list(filter(lambda x: (x != ''), output[1].strip().split(',')))

    if 'MCM_AAA' in domains_row:  # it's matrix_archaea_fast or matrix_bacteria_fast
        # deal with 'SAM_Arap_1,2,3_like'
        modified_domains_row = []
        for fname in domains_row:
            if fname.startswith('SAM_Arap'):
                modified_domains_row.append('SAM_Arap_1,2,3_like')
            elif fname == '2' or fname == '3_like':
                continue
            else:
                modified_domains_row.append(fname)

        # turn lists into dicts for faster lookup
        genomes_dict = {item: index for index, item in enumerate(genomes_row)}
        domains_dict = {item: index for index, item in enumerate(modified_domains_row)}

    else:
        # turn lists into dicts for faster lookup
        genomes_dict = {item: index for index, item in enumerate(genomes_row)}
        domains_dict = {item: index for index, item in enumerate(domains_row)}

    # parse the rest of output into list of lists
    output = output[2:]
    output = [list(map(int, filter(None, line.strip().split(',')))) for line in
              output]
    return domains_dict, genomes_dict, output


"""
3997 + ARUA distribution
"""

# file_path = 'xgroup_matrix_archaea_recovered.csv'
file_path = 'xgroup_matrix_bacteria_recovered.csv'
# file_path = 'xgroup_matrix_eukaryotes_recovered.csv'

domains_dict, genomes_dict, output = parseMatrix(file_path)
print(len(domains_dict))
print(len(genomes_dict))
#
# xgroups = ['3997', '1055', '136', '158', '185', '199', '2010', '239', '270', '3001', '3052', '309', '3115', '319', '3688', '374', '3754', '3896', '4004', '4019', '4022', '4028', '4048', '4076', '4110', '4126', '4160', '4237', '4295', '5084', '528', '582', '590', '6051', '6096', '640', '7513', '7551', '7553', '7567', '842', '869', '875', '876', '9']
# xgroups = ['3997', '1055', '158', '185', '199', '374', '4237', '582', '6051', '640', '869']  # excluded widely distributed xgroups

# xgroup2genomes = {}
# for x in xgroups:
#     for g in genomes_dict.keys():
#         if output[genomes_dict[g]][domains_dict[x]] > 0:  # xgroup exists in this genome
#             if x not in xgroup2genomes.keys():
#                 xgroup2genomes[x] = [g]
#             else:
#                 xgroup2genomes[x].append(g)
#     # record an empty list for xgroups that don't exist in any genome
#     if x not in xgroup2genomes.keys():
#         xgroup2genomes[x] = []
# print(len(xgroup2genomes))

# bacteria has too many genomes => let's just print the number of genomes for each xgroup
# xgroup2genomesNum = {}
# for key, lst in xgroup2genomes.items():
#     xgroup2genomesNum[key] = len(lst)

# dict2csv(xgroup2genomes, 'xgroup2genomes_archaea.csv')
# dict2csv(xgroup2genomes, 'xgroup2genomes_bacteria.csv')
# dict2csv(xgroup2genomes, 'xgroup2genomes_eukaryotes.csv')
# dict2csv(xgroup2genomesNum, 'xgroup2genomesNum_bacteria.csv')


"""
Josh annotation for aerobe, anaerobe, and facultative bacteria
"""

genome_path = 'gtdb_r207_3350genomesWithO2Requirements_copy.csv'
genome2oxy = csv2dict(genome_path)
print(len(genome2oxy))  # 3315

x2ns = csv2dict("data/xgroup2networkSize.csv")
print(len(x2ns))
xgroup397 = set(x2ns.keys())

xgroup2oxy = {}
for xgroup in xgroup397:  # 397
    xgroup2oxy[xgroup] = {'aerobe': 0, 'anaerobe': 0, 'facultative': 0, 'unknown': 0}
    for genome in genomes_dict.keys():  # ~62,000; contains 3-letter prefix
        if output[genomes_dict[genome]][domains_dict[xgroup]] > 0:  # xgroup exists in this genome
            if genome2oxy.get(genome[3:], "N/A") == ['Aerobe']:  # does not contain 3-letter prefix
                xgroup2oxy[xgroup]['aerobe'] += 1  # max 1711
            elif genome2oxy.get(genome[3:], "N/A") == ['Anaerobe']:
                xgroup2oxy[xgroup]['anaerobe'] += 1  # max 1187
            elif genome2oxy.get(genome[3:], "N/A") == ['Facultative']:
                xgroup2oxy[xgroup]['facultative'] += 1  # max 417

            else:
                xgroup2oxy[xgroup]['unknown'] += 1

print(xgroup2oxy['3997'])  # 42, 22, 15, 1112
print(xgroup2oxy['2003'])  # 1644, 1105, 409, 59133

dict2csv(xgroup2oxy, 'xgroup2oxy_bacteria.csv')


"""
genome2x instead
"""
# genome_path = 'gtdb_r207_3350genomesWithO2Requirements_copy.csv'
# genome2oxy = csv2dict(genome_path)
# print(len(genome2oxy))

# x2ns = csv2dict("data/xgroup2networkSize.csv")
# print(len(x2ns))
# xgroup397 = set(x2ns.keys())

# xgroup2230 = domains_dict.keys()

# genome2x = {}
# lost = set()
# for g in genome2oxy.keys():  # 3315 genomes
#     hits = {}
#
#     # take care of missing 3-letter prefix
#     pre_g = None
#     for genome in genomes_dict.keys():
#         if g in genome:
#             pre_g = genome  # found the genome
#             break
#
#     if pre_g is None:  # genome not found
#         lost.add(g)
#     else:
#         for x in xgroup2230:
#             hits[x] = output[genomes_dict[pre_g]][domains_dict[x]]
#         genome2x[g] = hits
#
# # print(len(genome2x))
# # print(len(lost))
# # dict2csv(genome2x, 'genome2x_bac_aerobeFull.csv')


"""
Eukaryotes hit distribution
"""
# file_path = 'xgroup_matrix_eukaryotes_recovered.csv'
#
# domains_dict, genomes_dict, output = parseMatrix(file_path)
# print(len(domains_dict))
# print(len(genomes_dict))
#
# xgroup2230 = domains_dict.keys()
#
# genome2x = {}
# for g in genomes_dict.keys():  # 196 genomes
#     hits = {}
#     for x in xgroup2230:
#         hits[x] = output[genomes_dict[g]][domains_dict[x]]
#     genome2x[g] = hits
#
# print(len(genome2x))
# dict2csv(genome2x, 'genome2x_euk_Full.csv')