#!/usr/bin/python -u

'''
Recurrence Rate is modeled as the derivative of the mean cumulative number of events through time.
The mean cumulative number of events through time is found from several sets (Monte Carlo) of 
potential cumulative number of events through time.
'''

import numpy as np
import matplotlib.pyplot as plt
import sys

### CHANGE THESE VARIABLES ###
minage,maxage,agestep = 0,350,10 #Ma
headerLines = 8 #For Straight-from-VEAM-Gui Output




###################
#Use all MC solutions to model Cumulative Number of Events
def cumulative_function(data,time):
	'''
	Creates a mean value of cumulative events for all times in the time list.
	Each row in the data array should be a list of event times. Each row is
	treated as one Monte Carlo solution of event times, for a Volcanic Event
	Age Model (VEAM)
	'''
	rows = len(data)                #returns (length,width) of data
	simulation_bins = range(rows+1) #histogram bins for the cumulative functions
	CDF_func = np.zeros(len(time))  #Mean Cumulative Function list
	
	for i,t in enumerate(time):
		#Time Histogram...
		#    count number of events before time t in each data row
		histogram_t = np.histogram(np.where(data>t)[0],bins=simulation_bins)[0]
		#[row 1: 0 events, row 2: 5 events, row 3: 2 events, ...]
	
		#Average of events in MC solutions with time (CDF)
		
		CDF_func[i] = np.mean(histogram_t)
	
	return CDF_func

###################
#Model RR by calculating derivative of mean Cumulative Distribution Function
def derivative(function):
	'''
	Calculates the derivative of a function using the Central difference method.
	Backward difference is used for the final element of the function,
	Forward difference is used for the first element.
	'''
	count = len(function)
	dfunction = np.zeros(count)
	for t in range(count):
		if t>0:
			if t<(count-1):
				dfunction[t] = (function[t-1]-function[t+1])/2 #Central  difference
			else:
				dfunction[t] = function[t-1]-function[t]       #Backward difference
		else:
			dfunction[t] = function[t]-function[t+1]         #Forward  difference
	return dfunction


###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################


#Get a list of events x 10k sets of potential ages
#MCdata = np.loadtxt("veam-results_arsia-mons_1klist.dat",skiprows=1,delimiter=",")

MCdata = np.loadtxt(sys.argv[1],skiprows=headerLines,delimiter=",")

(MCsets,eventCount) = np.shape(MCdata)
print "Loaded MC dataset"

times = np.arange(minage,maxage,agestep)

#################CUMULATIVE   UNCERTAINTY#######################
print "\n\n         ~~~~~~~~CDF of event count~~~~~~~~~"
#mean CDF of MC solutions
#CDF_all = np.load("CDF_each.npy")

CDF_each = np.zeros(( MCsets,len(times) ))
for i in range(len(MCdata)):
	CDF_each[i] = cumulative_function(MCdata[i:(i+1)],times)

#np.save("CDF_each",CDF_each)

#CDF_each = np.load("CDF_each5.npy")
print "Found cumulative event functions for each MC solution"


#################RECURRENCE RATE UNCERTAINTY#######################
print "\n\n         ~~~~~~~~Event RR~~~~~~~~~"

#"sawtooth" of each MC sol'n
RR_each = np.zeros(( MCsets,len(times) ))
for i,s in enumerate(CDF_each):
	RR_each[i] = derivative(s)/agestep

print "max RR:             ",np.max(RR_each)
print "min RR:             ",np.min(RR_each)
print "integrated mean RR: ",np.sum(np.mean(RR_each,axis=0))*agestep


#################CHARTS#####################
print "\n\n         ~~~~~~~~Charting~~~~~~~~~"
plt.figure(1)

xlimit = [minage,maxage]

###########CDF CHART#############
print "cumulative events"
plt.subplot(211)

simsPlot = min(1000,len(CDF_each))

for i in range(simsPlot):
	plt.plot(times,CDF_each[i],c='0.95')
plt.plot(times,CDF_each[0],c='0.95',label=('%d Individual Simulations' % simsPlot))
plt.plot(times,np.mean(CDF_each, axis=0),c='r',label='Mean Cumulative Events')
plt.plot(times,np.percentile(CDF_each, 10, axis=0),c='r',ls='--',label='10%')
plt.plot(times,np.percentile(CDF_each, 90, axis=0),c='r',ls='--',label='90%')


#Chart attributes
plt.xlim(xlimit)
plt.gca().invert_xaxis()
plt.ylabel('Cumulative Event Count')
plt.xlabel('Ma before present')
plt.title('Cumulative Events with time')
plt.legend(loc='upper left')


########RR CHART################
print "recurrence rate"
plt.subplot(212)
for i in range(simsPlot):
	plt.plot(times,RR_each[i],c='0.95')
plt.plot(times,RR_each[0],c='0.95',label=('%d Individual Simulations' % simsPlot))
plt.plot(times,np.mean(RR_each, axis=0),c='r',label="Mean Recurrence Rate")
plt.plot(times,np.percentile(RR_each, 10, axis=0),c='r',ls='--',label="10%")
plt.plot(times,np.percentile(RR_each, 90, axis=0),c='r',ls='--',label="90%")

#Chart attributes
plt.xlim(xlimit)
#plt.ylim([0,1.0])
plt.gca().invert_xaxis()
plt.ylabel('Events per Myr')
plt.xlabel('Ma before present')
plt.title('Recurrence Rate')
plt.legend(loc='upper left')


plt.show()



