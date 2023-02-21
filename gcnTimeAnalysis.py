import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams

mpl.rcParams.update({'font.size': 12, 'font.family': 'Times', 'mathtext.fontset': 'stix', 'lines.linewidth':2,'lines.markersize':5,\
                     'axes.linewidth':2, 'axes.labelsize': 'large', 'axes.labelweight': 'normal',\
                     'figure.figsize': (9.6, 6), 'ytick.major.size' : 6, 'xtick.major.size' : 6, \
                     'ytick.direction': 'in', 'xtick.direction': 'in', 'axes.labelpad': 15.0, \
                     'xtick.bottom': True, 'xtick.top': True, 'xtick.top': True, 'xtick.major.size': 9.5, \
                     'xtick.minor.bottom': False, 'xtick.minor.size': 5.0, 'xtick.minor.visible': True, \
                     'xtick.major.width': 1.8,'xtick.minor.width': 1.8,'ytick.left': True,'ytick.right': True, \
                     'ytick.minor.left': False, 'ytick.minor.right': False,'ytick.minor.size': 3.0, \
                     'ytick.major.width': 1.8,'ytick.minor.width': 1.8, 'ytick.minor.visible': True, \
                     'ytick.major.size': 9.5, 'axes.titlepad': 15.0})

import statistics as st
import numpy as np
import scipy as sy
import os
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit
from scipy import stats


with open('GRBGCNTimesDatabase.txt','r') as f:
    s = f.read()
    burstDict = eval(s)

#list of exceptional busts - 10 bursts with most circulars in the top 1000 delta_t
usualSuspects = ['130427A','211023A','100418A','060218','141121A','110328A','151027A','171205A','201015A','130831A']


times = [[] for i in range(2021-2006+1)]
timesNoRadio = [[] for i in range(2021-2006+1)]
timesNoSus = [[] for i in range(2021-2006+1)]
timesNoSusNoRadio = [[] for i in range(2021-2006+1)] #lists of delta_t segregated by year

lastTimes = [[] for i in range(2021-2006+1)]
lastTimesNoRadio = [[] for i in range(2021-2006+1)]
lastTimesNoSus = [[] for i in range(2021-2006+1)]
lastTimesNoSusNoRadio = [[] for i in range(2021-2006+1)] #lists of last observations (observations with highest delta_t)

minTimes = [[] for i in range(2021-2006+1)]
maxTimes = [[] for i in range(2021-2006+1)]  #list of lists of delta_t

numberOfBursts = [0 for i in range(2021-2006+1)]

#medianTime is a list of lists of gcn times since t_0. index 0 -> 2006; index 14 -> 2021

allTimes = [] #flat list of every t-t_0, when it's convenient to use
allTimesNoRadio = [] #as above, but with radio obs filtered out
allTimesNoSus = []
allTimesNoSusNoRadio = []
for (burstCode, gcnData) in burstDict.items(): #gcnData is a list of lists [[GCN, date time, delta_t, radioFlag]]
    burstYearIndex = int(burstCode[0:2])-6 #ex. '170817A' -> 17 -> index position 11
    numberOfBursts[burstYearIndex] += 1

    gcnData.sort(key=lambda gcn:gcn[2], reverse=True)

    lastFlag = False
    lastFlagNoRadio = False
    lastFlagNoSus = False
    lastFlagNoSusNoRadio = False
    for gcn in gcnData:
        times[burstYearIndex].append(gcn[2]) #list of t-t0 for every GCN, segregated by year
        allTimes.append(gcn[2]) #flat list of every t-t0 for every GCN
        if not lastFlag:
            lastTimes[burstYearIndex].append(gcn[2])
            lastFlag = True
        if not gcn[3]: #nonradio
            timesNoRadio[burstYearIndex].append(gcn[2]) #list of t-t0 for every GCN, segregated by year
            allTimesNoRadio.append(gcn[2]) #flat list of every t-t0 for every GCN
            if not lastFlagNoRadio:
                lastTimesNoRadio[burstYearIndex].append(gcn[2])
                lastFlagNoRadio = True
        if not (burstCode in usualSuspects): #no exceptional bursts
            timesNoSus[burstYearIndex].append(gcn[2]) #list of t-t0 for every GCN, segregated by year
            allTimesNoSus.append(gcn[2]) #flat list of every t-t0 for every GCN
            if not lastFlagNoSus:
                lastTimesNoSus[burstYearIndex].append(gcn[2])
                lastFlagNoSus = True
        if not (gcn[3] or (burstCode in usualSuspects)): #nonradio nonsuspect bursts
            timesNoSusNoRadio[burstYearIndex].append(gcn[2]) #list of t-t0 for every GCN, segregated by year
            allTimesNoSusNoRadio.append(gcn[2]) #flat list of every t-t0 for every GCN
            if not lastFlagNoSusNoRadio:
                lastTimesNoSusNoRadio[burstYearIndex].append(gcn[2])
                lastFlagNoSusNoRadio = True

for time in lastTimes:
    time.sort()
for time in lastTimesNoRadio:
    time.sort()
for time in lastTimesNoSus:
    time.sort()
for time in lastTimesNoSusNoRadio:
    time.sort()

# print(lastTimes[0])
# print(lastTimesNoRadio[0])
# print(lastTimesNoSus[0])
# print(lastTimesNoSusNoRadio[0])

print('Total of '+str(len(allTimes))+' times in the dataset')
print('Total of '+str(len(allTimesNoRadio))+' nonradio times in the dataset')
print('Total of '+str(len(allTimesNoSus))+' times in the nonexceptional dataset')
print('Total of '+str(len(allTimesNoSusNoRadio))+' nonradio times in the nonexceptional dataset')


#uncomment this block for the number of circulars per year
n = 2006
for year in times:  
    print('Total of '+str(len(year))+' times for '+str(n))
    n += 1

n = 2006
for year in timesNoRadio:  
    print('Total of '+str(len(year))+' nonradiotimes for '+str(n))
    n += 1

#uncomment this block for the number of bursts per year
for number in numberOfBursts:
  print(number)

medianTime = []
t10 = [] #time of first 10% of obs
t90 = [] #time if first 90% of obs
for time in times:
    time.sort()
    medianTime.append(st.median(time))
    t10.append(time[int(0.1*(len(time)-1))])
    t90.append(time[int(0.9*(len(time)-1))])


minimumTime = []
for time in times:
    minimumTime.append(min(time))

maximumTime = []
for time in times:
    maximumTime.append(max(time))

#plot of distribution of last observations


bulkTimes = [0 for year in times]
bulkTimesNoRadio = [0 for year in times]
bulkTimesNoSus = [0 for year in times]
bulkTimesNoSusNoRadio = [0 for year in times]



for yearIndex in range(len(lastTimes)):
    year = lastTimes[yearIndex]
    year.sort()
    flag = 0
    for index in range(len(year)):
        if (flag==0):
            fraction = index/len(year) #what fraction of bursts are we at?
            if fraction > 0.8: # -> 80th percentile of last observations
                bulkTimes[yearIndex] = year[index]
                flag += 1

for yearIndex in range(len(lastTimesNoRadio)):
    year = lastTimesNoRadio[yearIndex]
    year.sort()
    print(year)
    flag = 0
    for index in range(len(year)):
        if (flag==0):
            fraction = index/len(year) #what fraction of bursts are we at?
            if fraction > 0.8: # -> 80th percentile of last observations
                bulkTimesNoRadio[yearIndex] = year[index]
                flag += 1

for yearIndex in range(len(lastTimesNoSus)):
    year = lastTimesNoSus[yearIndex]
    year.sort()
    flag = 0
    for index in range(len(year)):
        if (flag==0):
            fraction = index/len(year) #what fraction of bursts are we at?
            if fraction > 0.8: # -> 80th percentile of last observations
                bulkTimesNoSus[yearIndex] = year[index]
                flag += 1

for yearIndex in range(len(lastTimesNoSusNoRadio)):
    year = lastTimesNoSusNoRadio[yearIndex]
    year.sort()
    flag = 0
    for index in range(len(year)):
        if (flag==0):
            fraction = index/len(year) #what fraction of bursts are we at?
            if fraction > 0.8: # -> 80th percentile of last observations
                bulkTimesNoSusNoRadio[yearIndex] = year[index]
                flag += 1

years = []
for i in range(2006,2021+1):
    years.append(i)

bulkTimeError = [0.1*time for time in bulkTimes] #10% error in the absense of a better measure for now
bulkTimeNoRadioError = [0.1*time for time in bulkTimesNoRadio]
bulkTimeNoSusError = [0.1*time for time in bulkTimesNoSus]
bulkTimeNoSusNoRadioError = [0.1*time for time in bulkTimesNoSusNoRadio]




#combined plots
plt.errorbar(years,[time/1000 for time in bulkTimes],yerr=[time/1000 for time in bulkTimeError],fmt="o",label='Complete dataset')
plt.errorbar(years,[time/1000 for time in bulkTimesNoRadio],yerr=[time/1000 for time in bulkTimeNoRadioError],fmt='o',label='Radio observations removed')
plt.xlabel('Year',fontsize=16)
plt.ylabel('$\Delta_t$ (ks)',fontsize=16)
plt.legend()
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
axes = plt.gca()
#axes.set_ylim([1,3.5e5])
plt.savefig('bulkLastTimeToPublish.eps', format='eps', dpi=1200)
plt.show()

#combined plot no suspects
plt.scatter(years,[time/1000 for time in bulkTimesNoSus],label='Complete dataset')
plt.scatter(years,[time/1000 for time in bulkTimesNoSusNoRadio],label='Radio observations removed')
plt.errorbar(years,[time/1000 for time in bulkTimesNoSus],yerr=[time/1000 for time in bulkTimeNoSusError],fmt="o")
plt.errorbar(years,[time/1000 for time in bulkTimesNoSusNoRadio],yerr=[time/1000 for time in bulkTimeNoSusNoRadioError],fmt='o')
plt.xlabel('Year',fontsize=16)
plt.ylabel('$\Delta_t$ (ks)',fontsize=16)
plt.legend()
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
axes = plt.gca()
#axes.set_ylim([1,3.5e5])
plt.savefig('bulkLastTimeToPublishNoSus.eps', format='eps', dpi=1200)
plt.show()

#linear-spaced bins plot of all t-t_0 for 2006-2021

allTimes.sort()
allTimes = allTimes[1:] #cropping out the one 0 time
allTimesNoRadio.sort()
allTimesNoRadio = allTimesNoRadio[1:] #cropping out the one 0 time

plt.hist([time for time in allTimes], bins = [i for i in range(0,100000000,1000)], histtype='step', lw=2,label='All times')
plt.hist([time for time in allTimesNoRadio], bins = [i for i in range(0,100000000,1000)], histtype='step', lw=2,label='No radio')
plt.xlabel('$\Delta_t$ (s)',fontsize=18)
plt.ylabel('Number of circulars',fontsize=18)
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.savefig('DistributionOfAllCircularsLin.eps',dpi=1200)

plt.show()

#log-spaced bins plot of all t-t_0 for 2006-2021

allTimes.sort()
allTimes = allTimes[1:] #cropping out the one 0 time
allTimesNoRadio.sort()
allTimesNoRadio = allTimesNoRadio[1:] #cropping out the one 0 time
plt.hist([np.log10(time) for time in allTimes], bins = int(len(allTimes)**0.5), histtype='step', lw=2,label='All times')
plt.hist([np.log10(time) for time in allTimesNoRadio], bins = int(len(allTimes)**0.5), histtype='step', lw=2,label='No radio')
plt.xlabel('Log($\Delta_t$) (s)',fontsize=18)
plt.ylabel('Number of circulars',fontsize=18)
plt.yscale('log')
plt.legend()
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.savefig('DistributionOfAllCircularsLog.eps',dpi=1200)

plt.show()

#plot of observations after t-t_0
allTimes.reverse()
allTimesNoRadio.reverse()
plt.plot([np.log10(time) for time in allTimes],[(index+1)/len(allTimes) for index in range(len(allTimes))],label='All times')
plt.plot([np.log10(time) for time in allTimesNoRadio],[(index+1)/len(allTimesNoRadio) for index in range(len(allTimesNoRadio))],label='No radio')

plt.xlabel('Log($\Delta_t$) (s)',fontsize=18)
plt.ylabel('F($\Delta_t$)',fontsize=18)
plt.grid(True,which='major',linewidth=1.5)
plt.grid(True,which='minor',ls='--')
plt.legend()
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.savefig('FractionAfterT-T0.eps',dpi=1200)

plt.show()


# #individual plots
# plt.errorbar(years,[time/1000 for time in bulkTimes],yerr=[time/1000 for time in bulkTimeError],fmt="o",label='80th percentile last observation')
# plt.xlabel('Year')
# plt.ylabel('$\Delta_t$ (ks)')
# plt.legend()
# axes = plt.gca()
# #axes.set_ylim([1,3.5e5])
# plt.show()

# plt.errorbar(years,[time/1000 for time in bulkTimesNoRadio],yerr=[time/1000 for time in bulkTimeNoRadioError],fmt="o",label='80th percentile last observation no radio')
# plt.xlabel('Year')
# plt.ylabel('$\Delta_t$ (ks)')
# plt.legend()
# axes = plt.gca()
# #axes.set_ylim([1,3.5e5])
# plt.show()

# plt.errorbar(years,[time/1000 for time in bulkTimesNoSus],yerr=[time/1000 for time in bulkTimeNoSusError],fmt="o",label='80th percentile last observation no usual suspects')
# plt.xlabel('Year')
# plt.ylabel('$\Delta_t$ (ks)')
# plt.legend()
# axes = plt.gca()
# #axes.set_ylim([1,3.5e5])
# plt.show()

# plt.errorbar(years,[time/1000 for time in bulkTimesNoSusNoRadio],yerr=[time/1000 for time in bulkTimeNoSusNoRadioError],fmt="o",label='80th percentile last observation no radio or usual suspects')
# plt.xlabel('Year')
# plt.ylabel('$\Delta_t$ (ks)')
# plt.legend()
# axes = plt.gca()
# #axes.set_ylim([1,3.5e5])
# plt.show()


# #plot of bulk time (t_80)
# bulkTimes = [0 for year in times]
# bulkTimeMinError = [0 for year in times]
# bulkTimeMaxError = [0 for year in times]

# for year in range(len(times)):
#     flag = 0
#     for index in range(len(times[year])):
#         if (flag==0):
#             fraction = index/len(times[year]) #what fraction of bursts are we at?
#             if fraction > 0.8: #0.8 -> t_80
#                 bulkTimes[year]=times[year][index]
#                 bulkTimeMinError[year]=bulkTimes[year]-times[year][index-int(index**0.5)] #poissonian uncertainty = sqrt of the count
#                 bulkTimeMaxError[year]=times[year][index+int(index**0.5)]-bulkTimes[year] #poissonian uncertainty = sqrt of the count
#                 flag += 1 #to break the loop

# years = []
# for i in range(2006,2021+1):
#     years.append(i)

# bulkTimeError = [bulkTimeMinError,bulkTimeMaxError]
# plt.scatter(years,bulkTimes)
# plt.errorbar(years,bulkTimes,yerr=bulkTimeError,fmt="o")
# plt.xlabel('Year')
# plt.ylabel('$\Delta_t$')
# axes = plt.gca()
# axes.set_ylim([1,3.5e5])
# plt.savefig('bulkTimeToPublish.eps', format='eps', dpi=1200)
# plt.show()



# #plot of distribution of circulars over time
# years = []
# for i in range(2006,2021):
#   years.append(i)

# timeErrorMin = []
# timeErrorMax = []

# for i in range(len(years)):
#   timeErrorMin.append(np.abs(medianTime[i]-t10[i]))
#   timeErrorMax.append(np.abs(medianTime[i]-t90[i]))

# timeError = [timeErrorMin,timeErrorMax]
# plt.scatter(years,medianTime)
# plt.errorbar(years,medianTime,yerr=timeError)
# plt.yscale("log")
# plt.xlabel('Year')
# plt.ylabel('Median time since t_0 (s)')
# axes = plt.gca()
# axes.set_ylim([1,3e6])
# plt.savefig('medianTimeToPublish.eps', format='eps', dpi=1200)
# plt.show()

# #plot of median GCN time vs. year

# plt.scatter(years,medianTime)
# plt.xlabel('Year')
# plt.ylabel('Median time since t_0 (s)')
# plt.title('Median time since t_0 vs. year')
# plt.show()

#slice up the data on a year-by-year basis and see if the distribution has changed?

# #first observations
# for year in years:
#   nBins=50
#   wBin=200 #width in seconds
#   bins=[]
#   for i in range(nBins):
#       bins.append(0) 

#   yearTimes = minTimes[int(year-2006)]
#   for time in yearTimes:
#       index = int(time//wBin)
#       if index >= nBins:
#           index=nBins-1
#       bins[index]+=1


#   xAxis = []
#   for i in range(nBins):
#       xAxis.append(i)

#   plt.scatter(xAxis,bins)
#   plt.xlabel('Time since t_0 (x200 s)')
#   plt.ylabel('GCN count')
#   axes = plt.gca()
#   axes.set_ylim([1,75])
#   plt.title(str(year)+' Earliest GCNs vs. time since t_0')
#   plt.show()

#last observations
# for year in years:
#   nBins=20
#   wBin=0.5 #width in log seconds
#   bins=[]
#   for i in range(nBins):
#       bins.append(0) 

#   yearTimes = maxTimes[int(year-2006)]
#   for time in yearTimes:
#       print(time)
#       index = int(time//wBin)

#       if index >= nBins:
#           index=nBins-1

#       bins[index]+=1


#   xAxis = []
#   for i in range(nBins):
#       xAxis.append(i)

#   plt.scatter(xAxis,bins)
#   plt.xlabel('Log time since t_0 (s) (x0.5)')
#   plt.ylabel('GCN count')
#   axes = plt.gca()
#   axes.set_ylim([1,40])
#   plt.title(str(year)+' Latest GCNs vs. time since t_0')
#   plt.show()

 # for year in years:
 #  nBins=200
 #  wBin=10000 #width in seconds
 #  bins=[]
 #  for i in range(nBins):
 #      bins.append(0) 
 #  yearTimes = times[int(year-2006)]
 #  for time in yearTimes:
 #      index = int(time//wBin)
 #      if index >= nBins:
 #          index=nBins-1
 #      bins[index]+=1


 #  xAxis = []
 #  for i in range(nBins):
 #      xAxis.append(i)
 #  plt.scatter(xAxis,bins)
 #  plt.xlabel('Time since t_0 (x10 ks)')
 #  plt.ylabel('GCN count')
 #  axes = plt.gca()
 #  axes.set_ylim([1,200])
 #  plt.title(str(year)+' GCNs vs. time since t_0')
 #  plt.show()

# #histogram of all observations since t_0

# for year in years:
#   nBins=100
#   wBin=10000 #width in seconds
#   bins=[]
#   for i in range(nBins):
#       bins.append(0) 

#   yearTimes = times[int(year-2006)]
#   for time in yearTimes:
#       index = int(time//wBin)

#       if index >= nBins:
#           index=nBins-1

#       bins[index]+=1


#   xAxis = []
#   for i in range(nBins):
#       xAxis.append(i)

#   plt.bar(xAxis,bins)
#   plt.xlabel('Time since t_0 (x10^5 s)')
#   plt.ylabel('GCN count')
#   axes = plt.gca()
#   #axes.set_ylim([1,40])
#   plt.title(str(year)+' All GCNs vs. time since t_0')
#   plt.show()

# #normalised histogram of all observations since t_0

# for year in years:
#   nBins=100
#   wBin=1e4 #width in seconds

#   histBins = []
#   for i in range(nBins):
#       histBins.append(i)

#   bins=[]
#   for i in range(nBins):
#       bins.append(0) 

#   yearTimes = times[int(year-2006)]
#   for time in yearTimes:
#       index = int(time//wBin)

#       if index >= nBins:
#           index=nBins-1

#       bins[index]+=1


#   xAxis = []
#   for i in range(nBins):
#       xAxis.append(i)

#   total = sum(bins)
#   normBins = bins
#   for i in range(len(normBins)):
#       normBins[i] = normBins[i]/total


#   if year==2016:
#       plt.hist([time/wBin for time in yearTimes],histBins,histtype='step', lw = 3)
#       plt.xlabel(r'Time since $t_0\,(\times 10^4\,\mathrm{s})$')
#       plt.ylabel('GCN count')
#       axes = plt.gca()
#       #axes.set_ylim([1,40])
#       #plt.title(str(year)+' normalised GCNs vs. time since t_0')
#       plt.savefig('2016NormalisedGCNSinceT0.eps',dpi=1200)
#       plt.show()

#histogram showing circulars vs. t-t_0 for sum of all years 2006-2021



# nBins=100
# wBin=1e4 #width in seconds

# histBins = [0 for i in range(nBins)]

# for time in allTimes:
#   index = int(time//wBin) #get the index of the observation

#   if index < nBins:
#       histBins[index] += 1

# buckets = []
# for i in range(nBins):
#   buckets.append(i)


# plt.hist(buckets,histBins,histtype='step', lw = 3)
# plt.xlabel(r'Time since $t_0\,(\times 10^4\,\mathrm{s})$')
# plt.ylabel('GCN count')
# axes = plt.gca()
# #plt.savefig('2016NormalisedGCNSinceT0.eps',dpi=1200)
# plt.show()

# #cumulative distribution
# for year in years:
#   nBins=1000
#   wBin=1000 #width in seconds
#   bins=[]
#   for i in range(nBins):
#       bins.append(0) 

#   yearTimes = times[int(year-2006)]
#   for time in yearTimes:
#       print(time)
#       index = int(time//wBin)

#       if index >= nBins:
#           index=nBins-1

#       bins[index]+=1

#   #producing cumulative distribution
#   cumulateBins = bins
#   for index in range(len(bins)):
#       cumulateBins[-(index+1)]=sum(bins[:-(index+1)])



#   timeSinceT0 = []
#   for i in range(nBins):
#       timeSinceT0.append(wBin*(i+1))

#   plt.scatter(timeSinceT0,cumulateBins)
#   plt.xlabel('Time since t_0 (ks)')
#   plt.ylabel('GCN count')
#   axes = plt.gca()
#   axes.set_ylim([1,1000])
#   plt.title(str(year)+' Cumulative GCNs vs. time since t_0')
#   plt.show()

# #normalised cumulative distribution with function
# def expCDF(x,lam):
#   return 1-np.exp(-x*lam)

# sampleYear1 = []
# sampleYear2 = []

# for year in years:
#   nBins=2000
#   wBin=1000 #width in seconds
#   bins=[]
#   for i in range(nBins):
#       bins.append(0) 

#   yearTimes = times[int(year-2006)]
#   for time in yearTimes:
#       index = int(time//wBin)

#       if index >= nBins:
#           index=nBins-1

#       bins[index]+=1

#   #producing normalised cumulative distribution
#   cumulateBins = [sum(bins) for bucket in bins]
#   total = sum(bins)
#   for index in range(len(bins)):
#       cumulateBins[index] = cumulateBins[index]-sum(bins[index+1:])
#       cumulateBins[index] = cumulateBins[index]/total


#   timeSinceT0 = []
#   for i in range(nBins):
#       timeSinceT0.append(wBin*(i+1))

#   p = plt.get_cmap('plasma')

#   if year%4==0:
#       plt.plot([x/1000 for x in timeSinceT0],cumulateBins,label=str(year),color=p((year-2006)/18))

#   if year == 2006:
#       sampleYear1 = cumulateBins
#   if year == 2018:
#       sampleYear2 = cumulateBins

#   # p0 = sy.array([0.000005])
#   # popt, pcov = curve_fit(expCDF, np.asarray(timeSinceT0), np.asarray(cumulateBins),p0)
#   # plt.plot([x/1000 for x in timeSinceT0],expCDF(np.asarray(timeSinceT0),*popt),label=str(year)+' model')




# plt.xlabel('Time since t_0 (ks)')
# plt.ylabel('GCN count')
# axes = plt.gca()
# axes.set_ylim([0,1])
# axes.set_xlim([0,1000])
# axes.legend()
# #plt.title('Normalised cumulative GCNs vs. time since t_0')
# plt.savefig('NormalisedCumulativeGCNs.eps', format='eps', dpi=1200)
# plt.show()


# print(str(stats.ks_2samp(sampleYear1,sampleYear2)))

#How long since t_0 until 80% of observations are made?

# bulkTime = [] #time until bulk of observations (80%) are made
# timeDev = [] #standard deviation of observation times for each year
# def expCDF(x,lam):
#   return 1-np.exp(-x*lam)



# for year in years:
#   nBins=2000
#   wBin=1000 #width in seconds
#   bins=[]
#   for i in range(nBins):
#       bins.append(0) 

#   yearTimes = times[int(year-2006)]
#   timeDev.append(st.stdev(yearTimes))

#   for time in yearTimes:
#       index = int(time//wBin)

#       if index >= nBins:
#           index=nBins-1

#       bins[index]+=1

#   #producing normalised cumulative distribution
#   cumulateBins = [sum(bins) for bucket in bins]
#   total = sum(bins)
#   for index in range(len(bins)):
#       cumulateBins[index] = cumulateBins[index]-sum(bins[index+1:])
#       cumulateBins[index] = cumulateBins[index]/total

#   flag = 0
#   for i in range(len(cumulateBins)):
#       if flag == 0 and cumulateBins[i] > 0.8: #0.8 -> t_80, etc.
#           bulkTime.append((i+1))
#           flag += 1


#   timeSinceT0 = []
#   for i in range(nBins):
#       timeSinceT0.append(wBin*(i+1))





# plt.scatter(years,bulkTime)
# plt.errorbar(years,bulkTime,yerr=np.sqrt(bulkTime),fmt='o')

# plt.xlabel('Year')
# plt.ylabel('Time since t_0 (ks)')
# axes = plt.gca()
# #axes.set_ylim([0, 300])
# # axes.set_xlim([200,1000])
# #plt.title('Time to publish 80 percent of GCNs')
# #plt.savefig('timeToPublish.eps', format='eps', dpi=1200)
# plt.show()


