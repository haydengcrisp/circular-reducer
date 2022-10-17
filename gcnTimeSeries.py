import matplotlib.pyplot as plt 
import statistics as st
import numpy as np
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


gcnDir = '/home/hayden/Desktop/paper2/unbiased GCN-based research/gcn3/'

startYear = 2006
endYear = 2020
timeEx = re.compile('([0-9]{2}:[0-9]{2}[:][0-9]{2})')
timeBodyEx = re.compile('([0-9]{2}:[0-9]{2}[:][0-9]{2}|[0-9]{2}:[0-9]{2}[:][0-9]{2}[.][0-9]+)(?:\s+\(*[Uu]|\(*[Uu])')
dateEx = re.compile('([0-9]{2}/[0-9]{2}/[0-9]{2})')
durationEx = re.compile('[0-9]+.[0-9]+( s| ks)')
t90Ex1 = re.compile('(?<=T90|t90).+') #returns all characters after T90 in a string
t90Ex2 = re.compile('[0-9]+( |\\.)*[0-9]*( s| ks| ms)') #matches things that look like times after T90

errors = []


for burstCode in allGRBs: #for each burst in my dataset
	t_0=''
	d_0=str(burstCode[0:2])+'/'+str(burstCode[2:4])+'/'+str(burstCode[4:6]) #date_0
	dt_0='' #datetime_0
	t_gcn=[] #list of lists, [GCN, date time, time since dt_0]
	firstGcnCode=1e10
	negativeTimeError=0
	noTriggerError=0
	noCircularError=0
	try:
		for gcn in os.listdir(gcnDir): #loops over every gcn in my archive to find all the relevant circulars
			with open(gcnDir+str(gcn),'r',encoding='latin-1') as f:
				circular = f.readlines() #list of lines in the circular
				if (len(circular) > 2) :#removes nonsense
					if (('GRB'+burstCode in circular[2]) or ('GRB '+burstCode in circular[2])) : #subject is the 3rd line in the header
	
	
							timePub=re.findall(timeEx,circular[3])
							datePub=re.findall(dateEx,circular[3]) #4th line in header is date/time
							dateString = datePub[0]+' '+timePub[0] #each should only contain one element
							t_gcn.append([gcn,dateString,0])
				f.close()
	
		print('All circulars for GRB '+str(burstCode)+' found!')

		if len(t_gcn)==0:
			noCircularError+=1
			raise(Exception)
			#sort by circular number to get chronological order of bursts
		t_gcn.sort(key=lambda tup: int(tup[0][:-5])) #strip terminating .gcn and convert to int before sorting



	#get the trigger time, which is almost always the only time in the body of the circular
		firstGcnCode = t_gcn[0][0]
		print('First circular: '+firstGcnCode)
		with open(gcnDir+str(firstGcnCode),'rb') as f:
			s = f.read()
			if (len(re.findall(timeBodyEx,s.decode('latin-1')))<1): #no trigger time?
				noTriggerError+=1
				f.close()
				raise(Exception)

			t_0 = re.findall(timeBodyEx,s.decode('latin-1'))[0] #1st time in report is the trigger
			print(burstCode+' triggered at UTC '+d_0+' '+t_0+', in '+str(firstGcnCode))
			dt_0 = dt.strptime(d_0+' '+t_0,'%y/%m/%d %H:%M:%S')
		f.close()
		
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

	except Exception:
		if (noCircularError > 0):
			errors.append('No circulars found related to GRB '+str(burstCode))

		if (noTriggerError > 0):
			errors.append('Index Error related to GCN '+str(firstGcnCode)+': no trigger time?')

		if (negativeTimeError>0):
			errors.append('Time error - '+gcn[0]+' published before trigger time of burst '+burstCode)


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

