# -*- coding: utf-8 -*-
import editdistance
import gzip
import sys

inFilename = sys.argv[1]
outFilename = sys.argv[2]

with gzip.open(inFilename, 'rt') as inFile: # handle gz
    with gzip.open(outFilename, 'w+t') as outFile:

            tripletDict = dict()   # triplet of gene, barcode combo (ie cell), UMI
            readDict = dict()
            for line in inFile:
                # revComp = None
                lineCols = line.split()
                gene = lineCols[1]
                if gene == '.':
                    continue
                # if lineCols[6] == 'NONE' and lineCols[8] == 'NONE':
                #     revComp = False
                #     cell = lineCols[1] + '-' + lineCols[3]
                #     umi = lineCols[5]
                # elif lineCols[1] == 'NONE' and lineCols[3] == 'NONE':
                #     revComp = True
                #     cell = lineCols[6] + '-' + lineCols[8]
                #     umi = lineCols[10]

                cell = lineCols[3]
                umi = lineCols[4]
                triplet = (gene, cell, umi)
                readDict[triplet] = lineCols[0]
                if triplet not in tripletDict:
                    tripletDict[triplet] = 1
                else:
                    tripletDict[triplet] += 1

            print('done counting triplets\nnumber of triplets:', len(tripletDict))
            # PROBABLY NEEDS TO BE REVERSE = TRUE?
            # below should be list of double tuples where first elem is triple tuple (triplet) and second elem. is count
            sortedTriplets = sorted(tripletDict.items(), key=lambda kv:(kv[1], kv[0]), reverse=True)
            # print(sortedTriplets[0], sortedTriplets[1], sortedTriplets[2])
            print('done ranking triplets')

            # make list consisting of (read, gene, cell, umi) preserving sorted order of sortedTriplets
            sortedTripletsWithRead = []

            for entry in sortedTriplets:
                sortedTripletsWithRead.append((readDict[entry[0]], entry[0][0], entry[0][1], entry[0][2]))

            # for entry in sortedTriplets:
                # tripletFile.write(readDict[entry[0]] + '\t' + '\t'.join(entry[0]) + '\n')
            # print('done writing triplet file')

            seenDoublets = dict()
            doubletCounter = 0
            notIn = 0
            yesIn = 0
            doubletCounterAlt = 0
            lineCount = 0

            approvedReads = {}

            # go through triplets and check for edit distances b/w UMIs w/ same doublet
            for trip in sortedTripletsWithRead:
                lineCount += 1 # meaning tripletCount

                ignoreUMI = False
                # doublet = (lineColsTrip[1], lineColsTrip[2])        # gene, cell
                doublet = (trip[1], trip[2])  # gene, cell

                doubletCounterAlt += 1
                if trip[3] == 'NONE':
                    continue
                if doublet not in seenDoublets:
                    notIn += 1
                    seenDoublets[doublet] = trip[3]
                    approvedReads[trip[0]] = 1
                    # doubletFile.write(lineColsTrip[0] + '\t' + lineColsTrip[1] + '\t' + lineColsTrip[2] + '\n')

                    doubletCounter += 1

                else:
                    yesIn += 1
                    # print('>1 UMI per gene–cell')

                    # choose method for edit_distance

                    for umi in seenDoublets[doublet].split():
                        if editdistance.eval(umi, trip[3]) < 4:
                            #print('small edit dist found')
                            ignoreUMI = True
                            break
                    if ignoreUMI:
                        continue

                    seenDoublets[doublet] += '\t' + trip[3]

                    # output format: read   gene    cell
                    approvedReads[trip[0]] = 1
                    # doubletFile.write(lineColsTrip[0] + '\t' + lineColsTrip[1] + '\t' + lineColsTrip[2] + '\n')

                    doubletCounter += 1


            inFile.seek(0)

            for inline in inFile:
                if inline.split()[0] in approvedReads:
                    outFile.write(inline)


            #print(seenDoublets)

            print('done writing doublet file\nnumber of doublets:', str(doubletCounter))
            print('notIn', notIn)
            print('yesIn', yesIn)
            print('doubletCounterAlt', doubletCounterAlt)
            print('lineCount', lineCount)
