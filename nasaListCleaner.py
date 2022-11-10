

with open('nasaListOfTransients.txt','r') as f:
	with open('nasaTempListOfGRBs.txt','w') as g:
		for entry in f.readlines():
			splitEntry = entry.split(' ')
			if splitEntry[0] == 'GRB':
				#minor weird terminators
				if '/' in splitEntry[1]:
					print(splitEntry[1][:splitEntry[1].rfind('/')],file=g)
				elif ':' in splitEntry[1]:
					print(splitEntry[1][:splitEntry[1].rfind(':')],file=g)
				#pre-2010 standardisation of multiple bursts on the same day
				elif ('A & B & C' in entry or 'A, B, and C' in entry or 'A and B and C' in entry):
					print(splitEntry[1]+'C',file=g)
					print(splitEntry[1]+'B',file=g)
					print(splitEntry[1],file=g)
				elif ('A & B' in entry or 'A and B' in entry):
					print(splitEntry[1][:6]+'B',file=g)
					print(splitEntry[1][:6],file=g)
				else:
					print(splitEntry[1],file=g)

	g.close()
f.close()

with open('nasaTempListOfGRBs.txt','r') as f:
	with open('nasaCleanListOfGRBs.txt','w') as g:
		for entry in f.readlines():
			if '(A,B)' in entry:
				print(entry[:6]+'B'.strip(),file=g)
				print(entry[:6].strip(),file=g)
			else:
				print(entry.strip(),file=g)
	g.close()
f.close()
