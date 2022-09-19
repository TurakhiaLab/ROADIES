with open('mapping.txt', 'r') as inFile:
    lines = inFile.read().split('\n')
    mapDict = dict()
    for line in lines[:-1]:
        [gene, spec] = line.split(' ')
        if spec in mapDict.keys():
            mapDict[spec].append(gene)
        else:
            mapDict[spec] = [gene]
    # print(len(mapDict.keys()))
    # for s in mapDict.keys():
    #     print(len(mapDict[s]))

with open('mapping2.txt', 'w') as outFile:
    for spec in mapDict.keys():
        line = spec + ': ' + ','.join(mapDict[spec]) + '\n'
        outFile.write(line)

    
