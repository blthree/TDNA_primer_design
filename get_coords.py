def lambda_handler(event, context):
    f = open('TAIR10_genes.tab', 'r')
    storage_dict = {}
    # build a dictionary of the start, end, orientation
    for line in f.readlines():
        s_line = line.strip('\n').split('\t')
        storage_dict[s_line[3]] = {'chr': s_line[0], 'start': int(s_line[1]), 'end': int(s_line[2]),
                                   'orientation': s_line[4]}

    # query dictionary for a specific gene
    try:
        result = storage_dict.get(event['gene_id'].upper())
    except KeyError:
        result = 'No match found'

    return result
    # raise Exception('Something went wrong')
