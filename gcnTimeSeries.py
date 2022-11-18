import matplotlib.pyplot as plt 
import statistics as st
import numpy as np
import pandas as pd
from datetime import datetime as dt
import os
import re
from scipy.interpolate import CubicSpline

with open('nasaCleanListOfGRBs.txt','r') as f:
	allGRBs = [entry.split('\n',maxsplit=1)[0] for entry in f.readlines()]
	f.close()

for i in range(len(allGRBs)):
	if (allGRBs[i][-1] == 'X' or allGRBs[i][-1] == 'S'): #X = no Gamma component, S = short; we don't care, as long as there's followup
		allGRBs[i]=allGRBs[i][:-1]

burstTimeDict = {}


gcnDir = '/home/hayden/Desktop/git-circular-reducer/circular-reducer/gcn3/'

startYear = 2006
endYear = 2021
timeEx = re.compile('([0-9]{2}:[0-9]{2}[:][0-9]{2})')
dateEx = re.compile('([0-9]{2}/[0-9]{2}/[0-9]{2})')

timeTrigEx1 = re.compile('([0-9]{1,2}:[0-9]{2}[:][0-9]{2}).*(?:\s+\(*[Uu]|\(*[Uu])') #matches 12:34:56 (UTC) and similar
timeTrigEx2 = re.compile('[0-9]{4}-[0-9]{2}-[0-9]{2}T([0-9]{2}:[0-9]{2}:[0-9]{2})') #matches 1999-12-31T12:34:56 and similar
timeTrigEx3 = re.compile('[0-9]{4}-[0-9]{2}-[0-9]{2}\s+([0-9]{2}:[0-9]{2}:[0-9]{2})') #matches 2021-12-26 05:25:41 and similar
timeTrigEx4 = re.compile('([0-9]{1,2}:[0-9]{2}:[0-9]{2})\s+[on]+\s+[0-9]+\s\w+') #matches 11:02:54 on 21 Oct and similar
timeTrigEx4_1 = re.compile('([0-9]{1,2}:[0-9]{2}:[0-9]{2}).[0-9]+\s[on]+\s\w+') #matches 11:02:54.00 on 21 Oct and similar
timeTrigEx5 = re.compile('[UTC]+\s\(*([0-9]{2}:[0-9]{2}:[0-9]{2}).*[0-9]*\)*') #matches UT (01:20:09) and similar
timeTrigEx6 = re.compile('([0-9]{1,2}:[0-9]{2})\s+[UTC]+') #matches partial times of the form 12:34 UTC and similar
timeTrigEx7 = re.compile('([0-9]{1,2})[hH]\s([0-9]{2})[mM]\s([0-9]{2})[sS]\s*[UTC]+') # matches 21h 13m 52s UTC and similar

errors = []


#circular cleaning block
#todo: stuff
def findPubTime(burst,circular): #burst is '160203A', circular is the string '1234.gcn3'
	with open(gcnDir+str(circular),'r',encoding='latin-1') as f:
		circularText = f.readlines()
		#check for faulty circulars
		if (len(circularText) < 4):
			f.close()
			return(None)
		else:
			baseCase = False
			masterCase = False
			preStandardCase = False
			
			subject = circularText[2].split()
			for word in subject:
				cleanWord = re.sub(r'[\W_]+', '', word).upper() #all caps, no special characters
				if (burst == cleanWord or 'GRB'+burst == cleanWord):
					baseCase = True
				elif (burst[:6]+'.' in word or 'GRB'+burst[:6]+'.' in word):
					masterCase = True
				elif (burst+'A'==cleanWord or (burst[-1]=='A' and burst[:6]==cleanWord)):
					preStandardCase = True

			if (baseCase or masterCase or preStandardCase) : #subject is the 3rd line in the header. header is 5 lines long
				timePub=re.findall(timeEx,circularText[3])
				datePub=re.findall(dateEx,circularText[3]) #4th line in header is date/time
				dateString = datePub[0]+' '+timePub[0] #each should only contain one element
				f.close()
				return(dateString)
			else:
				f.close()
				return(None)

# def findPubTime(burst,circular): #burst is '160203A', circular is the string '1234.gcn3'
# 	with open(gcnDir+str(circular),'r',encoding='latin-1') as f:
# 		circularText = f.readlines()
# 		if ((len(circularText) > 5) and (('GRB'+burst in circularText[2].upper()) or ('GRB '+burst in circularText[2].upper()))) : #subject is the 3rd line in the header. header is 5 lines long
# 			timePub=re.findall(timeEx,circularText[3])
# 			datePub=re.findall(dateEx,circularText[3]) #4th line in header is date/time
# 			dateString = datePub[0]+' '+timePub[0] #each should only contain one element
# 			f.close()
# 			return(dateString)
# 		else:
# 			f.close()
# 			return(None)

def findTriggerTime(circular):
	with open(gcnDir+circular,encoding='latin-1') as f:
		circularText = f.read()
		timeTrig=[]
		#we get the date from the name of the burst, so we only need to return the time in a standardised format
		#checking the time formats
		timeTrig=re.findall(timeTrigEx1,circularText) #12:34:56 (UTC) format
		if len(timeTrig)<1 : 
			timeTrig=re.findall(timeTrigEx2,circularText) #1999-12-31T12:34:56 format
		if len(timeTrig)<1 : 
			timeTrig=re.findall(timeTrigEx3,circularText) #2021-12-26 05:25:41 format
		if len(timeTrig)<1 : 
			timeTrig=re.findall(timeTrigEx4,circularText) #11:02:54 on 21 Oct format
		if len(timeTrig)<1 : 
			timeTrig=re.findall(timeTrigEx4_1,circularText) #11:02:54.00 on 21 Oct format
		if len(timeTrig)<1 : 
			timeTrig=re.findall(timeTrigEx5,circularText) #UT (01:20:09) format
		if len(timeTrig)<1 :
			timeTrigTemp = re.findall(timeTrigEx7,circularText) #21h 13m 52s UTC format, returns (21,13,52) and needs to be combined
			for time in timeTrigTemp:
				timeTrig.append(':'.join(time))
		#partial times
		if len(timeTrig)<1 : 
			tempTimeTrig=re.findall(timeTrigEx6,circularText) #12:34 UTC partial time; set seconds unit to 00
			for time in tempTimeTrig:
				timeTrig.append(time+':00')
		f.close()
		tempTimeTrig = timeTrig
		for time in tempTimeTrig:
			if time == '':
				timeTrig.remove(time)
	#special case: GRB 201014A trigger time taken from SWIFT trigger 1000255
	if burstCode == '201014A':
		timeTrig.append('22:48:38')
	#special case: GRB 200709B trigger time taken from HAWC notice 1009500
	if burstCode == '200709B':
		timeTrig.append('04:35:07')
	#special case: GRB 190404B trigger time taken from MAXI archive (http://maxi.riken.jp/grbs/)
	if burstCode == '190404B':
		timeTrig.append('13:14:29')
	#special case: GRB 181126B trigger time taken from FERMI trigger 564897175 
	if burstCode == '181126B':
		timeTrig.append('13:14:29')
	#special case: GRB 161224A trigger time taken from SWIFT trigger 728268 
	if burstCode == '161224A':
		timeTrig.append('22:13:34')
	#special case: GRB 140818A trigger time taken from GCN 16705
	if burstCode == '140818A':
		timeTrig.append('05:31:28')
	return(timeTrig[0])

def cleanListOfCirculars(listOfCirculars): #input: unsorted results of os.listdir(); outputs sorted list with non-.gcn3 files removed


	workingList = [filename for filename in listOfCirculars if '.gcn3' in filename] #don't iterate over the working list, also clears most of the error files
	for filename in listOfCirculars:
		if (not filename[:-5].isnumeric() and filename in workingList): #remove all the 'neg3.gcn3' and 'mistake.gcn3' circulars
			print('removed '+filename)
			workingList.remove(filename)


	workingList.sort(key=lambda circ: int(circ[:-5])) #strip terminating .gcn3 from circular and convert to int before sorting
	return(workingList)


allCircs = cleanListOfCirculars(os.listdir(gcnDir))


#dataset construction block


for burstCode in allGRBs: #for each burst in my dataset
	try:
		t_0=''
		d_0=str(burstCode[0:2])+'/'+str(burstCode[2:4])+'/'+str(burstCode[4:6]) #date_0
		dt_0='' #datetime_0
		t_gcn=[] #list of lists, [GCN, date time, time since dt_0]
		firstGcnCode=None
		negativeTimeError=0
		noTriggerError=0
		noCircularError=0
	
		print('Looking for circulars related to GRB '+burstCode)
		for circ in allCircs: #loops over every gcn in my archive to find all the relevant circulars
			if (findPubTime(burstCode,circ) != None):
				print((circ,findPubTime(burstCode,circ)))
				t_gcn.append([circ,findPubTime(burstCode,circ),0])
			#special cases
			#typo in subject line/unusual notation
			if (circ == '31021.gcn3' and burstCode == '211023B'): 
				print((circ,findPubTime('2121023B',circ)))
				t_gcn.append([circ,findPubTime('2121023B',circ),0])
			if (circ == '23957.gcn3' and burstCode == '190312A'): 
				print((circ,findPubTime('2121023B',circ)))
				t_gcn.append([circ,findPubTime('190312446',circ),0]) #BALROG notation


	
		print('All circulars for GRB '+str(burstCode)+' found!')
	
		#special cases
		if burstCode == '212102A' or burstCode == '180523B':
			print('Special case: '+burstCode+' does not exist')
			continue #this burst does not exist
	
		#sort by circular number to get chronological order of bursts
		t_gcn.sort(key=lambda tup: int(tup[0][:-5])) #strip terminating .gcn and convert to int before sorting
	
	
	
		#get the trigger time, which is almost always the only time in the body of the circular
		firstGcnCode = t_gcn[0][0]
		print('First circular: '+firstGcnCode)
		#special cases
		#followup circular published before detection circular
		if firstGcnCode == '30587.gcn3':
			firstGcnCode = '30589.gcn3'
		if firstGcnCode == '28690.gcn3':
			firstGcnCode = '28695.gcn3'
		if firstGcnCode == '24031.gcn3':
			firstGcnCode = '24041.gcn3'
		if firstGcnCode == '22882.gcn3':
			firstGcnCode = '22883.gcn3'
		if firstGcnCode == '22618.gcn3':
			firstGcnCode = '22619.gcn3'
		if firstGcnCode == '22155.gcn3':
			firstGcnCode = '22156.gcn3'
		if firstGcnCode == '20825.gcn3':
			firstGcnCode = '20826.gcn3'

		t_0 = findTriggerTime(firstGcnCode)
		print(d_0,t_0)
		dt_0 = dt.strptime(d_0+' '+t_0,'%y/%m/%d %H:%M:%S')
			
		for gcn in t_gcn: #loop over only the linked circulars
			with open(gcnDir+str(gcn[0]),'rb') as f:
				s = f.read() 
				times=re.findall(timeEx,s.decode('latin-1'))
				dates=re.findall(dateEx,s.decode('latin-1'))
				if (times!=[] and dates!=[]): #these GCNs are either mistakes, not observations, or ancient
					dateString = dates[0]+' '+times[0] #pulls the date and time from the header as a string
					deltaT = dt.strptime(dateString,'%y/%m/%d %H:%M:%S') - dt_0
					if deltaT.total_seconds() < 0: #check for nonsense
						negativeTimeError+=1
						f.close()
						raise(Exception)
						
					gcn[2] = deltaT.total_seconds()
				f.close()
	
		burstTimeDict[burstCode]=t_gcn

	except Exception as e:
			errors.append(e)
			if burstCode != None:
				errors.append(burstCode)
			if firstGcnCode != None:
				errors.append(firstGcnCode)
			if negativeTimeError>0:
				errors.append(gcn[0]+'-ve time error')



#the practice of ending all bursts with a letter didn't start until ~2009, so we need to prune double entries
#i.e., removing circulars from GRB yymmddB from the entry for GRB yymmdd

for burst, gcnData in burstTimeDict.items(): #over all bursts
	for otherBurst, otherGcnData in burstTimeDict.items(): #check every other burst:
		if (burst != otherBurst and burst in otherBurst): #burst is the first event (GRB yymmdd), otherBurst is the later events (GRB yymmddB)
			for circular in gcnData:
				for otherCircular in otherGcnData: #checking every circular with every other circular
					if circular[0]==otherCircular[0]: #entry is in both
						gcnData.remove(circular) #remove the late burst's circular from the first burst's list
						print('gcn '+str(otherCircular[0])+' was in '+burst+' and '+otherBurst)



#write the whole thing to a file so i don't need to run this O(n^2) search every time i want a graph
with open('GRBGCNTimesDatabase.txt','w') as data:
	data.write(str(burstTimeDict))

with open('GRBGCNTimesDatabaseErrors.txt','w') as data:
	data.write(str(errors))