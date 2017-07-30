import pymongo
def test():
    connection_str = "mongodb://blthree:4o#a2k28@cluster0-shard-00-00-mbk9i.mongodb.net:27017,cluster0-shard-00-01-mbk9i.mongodb.net:27017,cluster0-shard-00-02-mbk9i.mongodb.net:27017/testDB?ssl=true&replicaSet=Cluster0-shard-0&authSource=admin"

    client = pymongo.MongoClient(connection_str)
    db = client.testDB.salk_table
    a = db.find({'name': 'SALK_105079'})

    str_docs = []
    for doc in a:
        print(doc)



def load_into_mongo():
    records = {}
    f = open('T-DNA.SALK', 'r')
    for line in f.readlines():
        s_line = line.strip('\n').split('\t')
        stock_name = s_line[0].split('.')[0]
        poly_name = s_line[0]
        poly_chr = s_line[1].split(':')[0]
        poly_chr = poly_chr[3:]
        poly_locs = s_line[3].split(',')[0]
        poly_start = poly_locs.split('/')[1].split('-')[0]
        poly_end = poly_locs.split('/')[1].split('-')[1]
        orientation = s_line[3].split('/')[0]
        # escape . in poly_name by replacing with ^
        poly_name = poly_name.replace('.', '^')
        # load into dict of dicts object
        if stock_name in records and poly_name in records[stock_name]['T-DNA insertions']:
            records[stock_name]['T-DNA insertions'][poly_name + '^1'] = {
                'Coordinates': {'chr': poly_chr, 'orientation': orientation,
                                'start': poly_start,
                                'end': poly_end}}
        elif stock_name not in records:
            records[stock_name] = {'T-DNA insertions': {poly_name: {
                'Coordinates': {'chr': poly_chr, 'orientation': orientation, 'start': poly_start, 'end': poly_end}}}}
        else:
            records[stock_name]['T-DNA insertions'][poly_name] = {
                'Coordinates': {'chr': poly_chr, 'orientation': orientation, 'start': poly_start,
                                'end': poly_end}}
    a = []
    for key in records:
        b = records[key]
        b['name'] = key
        a.append(b)
    db.insert_many(a)



connection_str = "mongodb://blthree:4o#a2k28@cluster0-shard-00-00-mbk9i.mongodb.net:27017,cluster0-shard-00-01-mbk9i.mongodb.net:27017,cluster0-shard-00-02-mbk9i.mongodb.net:27017/testDB?ssl=true&replicaSet=Cluster0-shard-0&authSource=admin"
client = pymongo.MongoClient(connection_str)
db = client.salkDB.salk_table
a = db.find({'name': 'SALK_000030'})
str_docs = []
print(a[0])
j = a[0]['T-DNA insertions']
print(j)
