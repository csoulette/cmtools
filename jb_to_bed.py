import os, sys

with open(sys.argv[1]) as lines:
    header = next(lines).rstrip().split("\t")
    for i in lines:
        cols = i.rstrip().split("\t")

        asType = cols[1]

        if asType == "cassette":
            inclJuncs = cols[6]
            data = inclJuncs.split(";")
            leftJuncs, rightJuncs = data[0].split(",")[0],data[1].split(",")[0]
            strand, gene = cols[4], cols[2]
            
            chrom,coords = leftJuncs.split(":")
            lc1, lc2 = [int(x) for x in coords.split("-")]

            chrom,coords = rightJuncs.split(":")
            rc1, rc2 = [int(x) for x in coords.split("-")]

            coords = sorted([rc1,rc2,lc1,lc2])

            dpsi = float(cols[-3])
            adjp = float(cols[-1])
            if adjp < 0.05:
                if dpsi>=30:
                    color = "165,0,38"
                elif dpsi>=10:
                    color = "244,109,67"
                
                else:
                    color = "254,224,144"
            else:
                if dpsi>10:
                    color = "224,243,248"
                else:
                    color = "189,189,189"
            sizes = "50,%s,50" % (abs(coords[2]-coords[1])-1)
            starts = "0,%s,%s" % ((coords[1]-coords[0])+50,(coords[-1]-coords[-2])+abs(coords[2]-coords[1])+(coords[1]-coords[0])+50)
            print(chrom, coords[0]-50, coords[-1]+50, "%s;%s" % (cols[8],gene), 0, strand, coords[0]-50, coords[-1]+50, color, 3, sizes, starts, sep="\t")
