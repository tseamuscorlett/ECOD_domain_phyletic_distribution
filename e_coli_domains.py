import csv

def dict2csv(dict, csv_file_path):
    # Write the dictionary to a CSV file
    with open(csv_file_path, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)

        # Write the header
        csv_writer.writerow(dict.keys())

        # Write the data
        try:
            csv_writer.writerows(zip(*dict.values()))
        except:
            csv_writer.writerows(dict.values())


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


ecoli_genome = 'RS_GCF_003697165.2'
arc_genome = 'GB_GCA_014729345.1'

file_path = 'xgroup_matrix_bacteria_recovered.csv'
# file_path = 'matrix_bacteria_recovered.csv'
# file_path = 'matrix_archaea_fast.csv'

domains_dict, genomes_dict, output = parseMatrix(file_path)

print(len(domains_dict))
g = genomes_dict[ecoli_genome]
xgroup2freq_ecoli = {}
for d in domains_dict.keys():
    # print(d, domains_dict[d], output[g][domains_dict[d]])
    xgroup2freq_ecoli[d] = [output[g][domains_dict[d]]]

print(xgroup2freq_ecoli)

dict2csv(xgroup2freq_ecoli, 'xgroup2freq_ecoli.csv')


