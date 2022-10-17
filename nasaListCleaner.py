with open('nasaListOfTransients.txt','r') as f:
	with open('nasaCleanListOfGRBs.txt','w') as g:
		for entry in f.readlines():
			splitEntry = entry.split(' ')
			print(splitEntry[0])
			if splitEntry[0] == 'GRB':
				print(splitEntry[1][:-1],file=g)

	g.close()
f.close()