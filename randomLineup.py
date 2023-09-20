import numpy as np
import random

with open('GRBGCNTimesDatabase.txt','r') as f:
    s = f.read()
    burstDict = eval(s)
    f.close()

randomBursts = [] #list of key-value pairs
for burst in range(20): #seems like a reasonable number
    randomBursts.append(random.choice(list(burstDict)))


with open('randomBursts.txt','w') as data:
    data.write(str(randomBursts))
    