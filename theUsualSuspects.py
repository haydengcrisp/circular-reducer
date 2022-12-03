import numpy as np

with open('GRBGCNTimesDatabase.txt','r') as f:
    s = f.read()
    burstDict = eval(s)

suspectList = [] #each entry is a list of the form [GRB, circular, log10(delta_t)]
for key in burstDict:
    for circ in burstDict[key]:
        suspect = [key, circ[0],np.log10(circ[2])]
        suspectList.append(suspect)

suspectList.sort(key=lambda suspect: suspect[2])
suspectList.reverse()

i = 0
while i < 100:
    print(suspectList[i])
    i += 1