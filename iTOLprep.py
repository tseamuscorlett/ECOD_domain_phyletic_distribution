import csv

def dict2csv(dict, csv_file_path):
    # Write the dictionary to a CSV file
    with open(csv_file_path, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)

        # Write the header
        csv_writer.writerow(dict.keys())

        # Write the data
        csv_writer.writerow(dict.values())


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


# file_path = 'xgroup_matrix_archaea_recovered.csv'
file_path = 'xgroup_matrix_bacteria_recovered.csv'
# file_path = 'xgroup_matrix_eukaryotes_recovered.csv'


domains_dict, genomes_dict, output = parseMatrix(file_path)
print(len(domains_dict))
#
xgroups = ['3997', '1055', '136', '158', '185', '199', '2010', '239', '270', '3001', '3052', '309', '3115', '319', '3688', '374', '3754', '3896', '4004', '4019', '4022', '4028', '4048', '4076', '4110', '4126', '4160', '4237', '4295', '5084', '528', '582', '590', '6051', '6096', '640', '7513', '7551', '7553', '7567', '842', '869', '875', '876', '9']
xgroups = ['3997', '1055', '158', '185', '199', '374', '4237', '582', '6051', '640', '869']

xgroup2genomes = {}
for x in xgroups:
    for g in genomes_dict.keys():
        if output[genomes_dict[g]][domains_dict[x]] > 0:  # xgroup exists in this genome
            if x not in xgroup2genomes.keys():
                xgroup2genomes[x] = [g]
            else:
                xgroup2genomes[x].append(g)
    # record an empty list for xgroups that don't exist in any genome
    if x not in xgroup2genomes.keys():
        xgroup2genomes[x] = []
print(len(xgroup2genomes))


# bacteria has too many genomes => let's just print the number of genomes for each xgroup
# xgroup2genomesNum = {}
# for key, lst in xgroup2genomes.items():
#     xgroup2genomesNum[key] = len(lst)


# dict2csv(xgroup2genomes, 'xgroup2genomes_archaea.csv')
dict2csv(xgroup2genomes, 'xgroup2genomes_bacteria.csv')
# dict2csv(xgroup2genomesNum, 'xgroup2genomesNum_bacteria.csv')
# dict2csv(xgroup2genomes, 'xgroup2genomes_eukaryotes.csv')


