#!/usr/bin/python -u

'''
Recurrence Rate is modeled as the derivative of the mean cumulative number of events through time.
The mean cumulative number of events through time is found from several sets (Monte Carlo) of 
potential cumulative number of events through time.
'''

import numpy as np
import matplotlib.pyplot as plt
import sys

minage,maxage = 0,350

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
#Use all MC solutions to model Cumulative Number of Volumes
def cumulative_volume(data,time,volumes):
	'''
	Creates a mean value of volume emitted before each time in the time list.
	Each row in the data array should be a list of event times. Each row is
	treated as one Monte Carlo solution of event times, for a Volcanic Event
	Age Model (VEAM)
	List of volumes must be ordered corresponding to the data array columns
	'''
	rows = len(data)                #returns (length,width) of data
	simulation_bins = range(rows+1) #histogram bins for the cumulative functions
	vol_cdf_func = np.zeros(len(time))  #Mean Cumulative Function list
	
	#if multiple MC soln's are being used at once to calculate vol flux
	if len(np.shape(data))>1:
		for i,t in enumerate(time):
			vol_cdf_func[i] = np.sum(volumes[np.where(data>t)[1]])/rows #gives mean volume output per MC set
	
	#if only one age set is being used for volume flux calculation
	else:
		for i,t in enumerate(time):
			vol_cdf_func[i] = np.sum(volumes[np.where(data>t)[0]]) #gives mean volume output per MC set
			
	return vol_cdf_func

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

MCdata = np.loadtxt(sys.argv[1],skiprows=1,delimiter=",")

(MCsets,eventCount) = np.shape(MCdata)
print "Loaded MC dataset"

tstp = 10
tmax = 500
tmin = 0
time5 = np.arange(tmin,tmax,tstp)

#################CUMULATIVE   UNCERTAINTY#######################
print "\n\n         ~~~~~~~~CDF of event count~~~~~~~~~"
#mean CDF of MC solutions
#CDF_all = np.load("CDF_each.npy")

CDF_each = np.zeros(( MCsets,len(time5) ))
for i in range(len(MCdata)):
	CDF_each[i] = cumulative_function(MCdata[i:(i+1)],time5)

#np.save("CDF_each",CDF_each)

#CDF_each = np.load("CDF_each5.npy")
print "Found cumulative event functions for each MC solution"


#################RECURRENCE RATE UNCERTAINTY#######################
print "\n\n         ~~~~~~~~Event RR~~~~~~~~~"

#"sawtooth" of each MC sol'n
RR_each = np.zeros(( MCsets,len(time5) ))
for i,s in enumerate(CDF_each):
	RR_each[i] = derivative(s)/tstp

print "max RR:             ",np.max(RR_each)
print "min RR:             ",np.min(RR_each)
print "integrated mean RR: ",np.sum(np.mean(RR_each,axis=0))*tstp



#################CUMULATIVE VOLUME UNCERTAINTY#####################
print "\n\n         ~~~~~~~~Incorporating Volumes~~~~~~~~~"
volumes = np.loadtxt("edifice-area_sorted-as-veam.dat",skiprows=1)
flow_thickness = 0.01 # km
volumes *= flow_thickness

print "Flow thickness assumed to be ",flow_thickness," km."

cum_vols = np.zeros(( MCsets,len(time5) ))
for i in range(MCsets):
	cum_vols[i] = cumulative_volume(MCdata[i],time5,volumes)
	
print "Total volume: ",np.max(cum_vols),np.sum(volumes), "(should be same no.)"


#################VOLUME FLUX UNCERTAINTY#####################
print "\n\n         ~~~~~~~~Volume Flux over time!~~~~~~~~~"

#cum_vols
#np.save("CDF_each",CDF_each)

#CDF_each = np.load("CDF_each5.npy")
#print "Found cumulative event functions for each MC solution"

#"sawtooth" of each MC sol'n
VF_each = np.zeros(( MCsets,len(time5) ))
for i,s in enumerate(cum_vols):
	VF_each[i] = derivative(s)/tstp

print "max VF:             ",np.max(VF_each)
print "min VF:             ",np.min(VF_each)
print "integrated mean VF: ",np.sum(np.mean(VF_each,axis=0))*tstp


#################CHARTS#####################
print "\n\n         ~~~~~~~~Charting~~~~~~~~~"
plt.figure(1)

xlimit = [minage,maxage]

###########CDF CHART#############
print "cumulative events"
plt.subplot(221)

for i in range(1000):
	plt.plot(time5,CDF_each[i],c='0.95')
plt.plot(time5,np.mean(CDF_each, axis=0),c='r',label='Mean Cumulative Events')
plt.plot(time5,np.percentile(CDF_each, 10, axis=0),c='r',ls='--',label='10%')
plt.plot(time5,np.percentile(CDF_each, 90, axis=0),c='r',ls='--',label='90%')


#Chart attributes
plt.xlim(xlimit)
plt.gca().invert_xaxis()
plt.ylabel('Cumulative Event Count')
plt.xlabel('Ma before present')
plt.title('Cumulative Events with time')
plt.legend(loc='upper left')


########RR CHART################
print "recurrence rate"
plt.subplot(222)
for i in range(1000):
	plt.plot(time5,RR_each[i],c='0.95')
plt.plot(time5,np.mean(RR_each, axis=0),c='r',label="Mean Recurrence Rate")
plt.plot(time5,np.percentile(RR_each, 10, axis=0),c='r',ls='--',label="10%")
plt.plot(time5,np.percentile(RR_each, 90, axis=0),c='r',ls='--',label="90%")

#Chart attributes
plt.xlim(xlimit)
plt.ylim([0,1.0])
plt.gca().invert_xaxis()
plt.ylabel('Events per Myr')
plt.xlabel('Ma before present')
plt.title('Recurrence Rate')
plt.legend(loc='upper left')


###########CUM_V CHART###########
print "cumulative volume"
plt.subplot(223)

for i in range(1000):
	plt.plot(time5,cum_vols[i],c='0.95')
plt.plot(time5,np.mean(cum_vols, axis=0),c='r',label='Mean Cumulative Volume')
plt.plot(time5,np.percentile(cum_vols, 10, axis=0),c='r',ls='--',label='10%')
plt.plot(time5,np.percentile(cum_vols, 90, axis=0),c='r',ls='--',label='90%')

#Chart attributes
plt.xlim(xlimit)
plt.gca().invert_xaxis()
plt.ylabel('Cumulative Volume, km^3')
plt.xlabel('Ma before present')
plt.title('Total volume erupted, with %0.2f km thick flows' % flow_thickness)
plt.legend(loc='upper left')


###########VF CHART###########
print "volume flux"
plt.subplot(224)

for i in range(1000):
	plt.plot(time5,VF_each[i],c='0.95')
plt.plot(time5,np.mean(VF_each, axis=0),c='r',label='Mean Volume Flux')
plt.plot(time5,np.percentile(VF_each, 10, axis=0),c='r',ls='--',label="10%")
plt.plot(time5,np.percentile(VF_each, 90, axis=0),c='r',ls='--',label="90%")

#Chart attributes
plt.xlim(xlimit)
plt.gca().invert_xaxis()
plt.ylabel('km^3 per Myr')
plt.xlabel('Ma before present')
plt.title('Volume Flux, with %0.2f km thick flows' % flow_thickness)
plt.legend(loc='upper left')



plt.show()



