from ssCorrect import BED12
import sys, os
from tqdm import *


first = BED12(sys.argv[1])
a = set()
for line in tqdm(first.getLine()):
	juncs = first.bed12toJuncs()
	if len(juncs)<1: continue
	a.add(tuple(juncs))

#print(len(a))

second = BED12(sys.argv[2])
b = set()
c = set()
l = dict()
for line in tqdm(second.getLine()):
	juncs = second.bed12toJuncs()
	if len(juncs)<1: continue
	juncs = tuple(juncs)
	b.add(juncs)
	l[juncs] = second.name#second.name.split("_")[0]
	if "ENST" in second.name:
		c.add(juncs)

#print(len(b))
#print("sig overlap",len(a.intersection(b)))
##print("sig annot",len(b.intersection(c)))
print("short annot", len(a.intersection(c)))

#[print(x) for x in list(l.values())]

[print(l[x]) for x in list(b.difference(a)) if "UPP1" in l[x]]

print("excl")

[print(x) for x in list(l.values()) if "UPP1" in x]

print("incl")

[print(l[x]) for x in list(b.intersection(a)) if "UPP1" in l[x]]
