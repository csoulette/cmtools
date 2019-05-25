import os, sys

with open(sys.argv[1]) as lines:
    for line in lines:
        
        c = line.rstrip().split()
        
        
        
        
        if c[4] == "+":
            c[1] = str(int(c[1])-50)
            d = int(c[2]) - int(c[1])
            print("\t".join(c), c[1], c[2], "0,0,0", "2", "50,50,", "0,%s," % (d), sep="\t")
        else:
            d = int(c[2]) - int(c[1]) + 50
            c[1] = str(int(c[1])-50)
            c[2] = str(int(c[2])+50)

            print("\t".join(c), c[1], c[2], "0,0,0", "2", "50,50,", "0,%s," % (d), sep="\t")
