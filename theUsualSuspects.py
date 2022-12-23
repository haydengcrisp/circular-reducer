import numpy as np

with open('GRBGCNTimesDatabase.txt','r') as f:
    s = f.read()
    burstDict = eval(s)

suspectList = [] #each entry is a list of the form [GRB, offenceCount, circular, log10(delta_t)]
for key in burstDict:
    for circ in burstDict[key]:
        suspect = [key, circ[0],np.log10(circ[2])]
        suspectList.append(suspect)

suspectList.sort(key=lambda suspect: suspect[2])
suspectList.reverse()

topSuspects = suspectList[:1000]
for suspect in topSuspects:
    count = 0
    burst = suspect[0]
    for otherSuspect in topSuspects:
        if burst == otherSuspect[0]:
            count += 1
    suspect.append(count)


i = 0
print('Burst, offending circular, log_10(delta_t), number of circs in top 1000')
while i < 1000:
    print(topSuspects[i])
    i += 1

topOffenders = []
for suspect in topSuspects:
    if suspect[3] > 1 and not suspect[0] in topOffenders:
        topOffenders.append(suspect[0])

print('Ordered by number of circulars in top 1000')
topSuspects.sort(key=lambda suspect: suspect[3])
topSuspects.reverse()
for suspect in topSuspects:
    print(suspect)
