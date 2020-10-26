#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 16:42:04 2017
@author: James Wilson, Jacob Richardson

VEAM is a VOLCANIC EVENT AGE MODELER
"""

import sys, os, time, operator
import threading, queue
from PyQt5 import QtWidgets, QtCore, QtGui
import numpy as np
from scipy.stats import truncnorm
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

global UINT_MAX
UINT_MAX = np.iinfo('uint32').max

### VEAM FUNCTIONS ###

def load_databases(strat_db_file, ages_db_file):
	#load age database using genfromtxt into 2xN string matrix
	try:
		relationships = np.genfromtxt(strat_db_file,skip_header=1,delimiter=',',dtype='unicode')
	except ValueError:
		sys.stderr.write('\n ERROR: Check for extra spaces, commas in the relationships in %s and try again\n' % ages_db_file)
		return -1, None
	#load ages database using genfromtxt
	try:
		AgeUncertainty = np.genfromtxt(ages_db_file, skip_header=1,delimiter=',', dtype='unicode')
	except:
		ValueError
		sys.stderr.write('\n ERROR: Check for extra spaces, commas in the names of events in %s and try again\n' % ages_db_file)
		return -1, None
	'''
	sys.stdout.write('\nThis is the age database\n', AgeUncertainty)
	'''
	
	eventLib = eventLibrary()
	
	# Add events to field
	for i,entry in enumerate(AgeUncertainty):
		exist = 0
		idname = entry[0]
		age    = float(entry[1])
		uncert = float(entry[2])
		
		if len(eventLib.events) > 0:
			for e in eventLib.events:
				if e.id==idname: #if event already exists, just add the new age model
					e.addAgeModel(age, uncert)
					exist = 1
					break
		if exist == 0:
			eventLib.addEvent(idname)
			eventLib.events[-1].addAgeModel(age, uncert)
	
	sys.stdout.write('      Loaded all Events from Ages Database successfully!\n')
	
	
	
	#Add Stratigraphic Relationships to Events
	for link in relationships:
		#The lower event is link[0], the higher event is link[1]
		#event.stratAbove and stratBelow
		both = 0
		for event in eventLib.events:
			#if the event is the lower event, append the higher event to stratAbove
			if event.id == link[0]:
				event.stratAbove.append(link[1])
				both += 2
			#if the event is the higher event, append the higher event to stratBelow
			elif event.id == link[1]:
				event.stratBelow.append(link[0])
				both += 1
			#If both events have been found, go to next relationship
			if both == 3:
				break
		if both == 0:#if both still = 0, neither event found
			sys.stderr.write('Events \'%s\' and \'%s\' in stratigraphic database not found in Age database!\n' %
																				(link[0],link[1]))
			return -1, None
		elif both != 3: #if both still = 1, lower event not found, if = 2, higher event not found
			sys.stderr.write('Event \'%s\' in Stratigraphic database not found in Age database!\n' % 
																				(link[both-1]))
			return -1, None
	
	
	#Check Stratigraphy then find expanded stratigraphic relationships for all events
	check = eventLib.checkAllStrat(relationships)
	if check == -1:
			return -1, None
	check = eventLib.getFullStratLists(len(relationships))
	if check == -1:
			return -1, None
	sys.stdout.write('Loaded all Relationships in Stratigraphy Database successfully!\n')
	
	return 0, eventLib

def sample_ages(eventsLib):
	'''
	This function generates a dictionary of sampled ages. 
	It iterates over a list of events. 
	It finds an acceptable age range based upon stratigraphic relationships and SampledAges. 
	This acceptable age range then truncates the probability distribution function. 
	In cases where minage/maxage from stratigraphy truncate the probability distribution function 
		to one of the tails where the probability of an event is zero, the code then chooses a 
		random uniform age between the two bounding units. 
	In this way, stratigraphy can override the radiometric date.
	'''
	
	### Give each event a new age and initialize it to be -1
	for event in eventsLib.events:
		event.veamAges.append(-9999)
		
	for event in eventsLib.events:
		#####FIND EVENT AGE RANGE######
		AcceptableAge_MIN = GLOBAL_MINAGE
		if len(event.allAboveInd) > 0:
			for e in event.allAboveInd:
				if eventsLib.events[e].veamAges[-1] > AcceptableAge_MIN:
					AcceptableAge_MIN = eventsLib.events[e].veamAges[-1]
		
		AcceptableAge_MAX = GLOBAL_MAXAGE
		if len(event.allBelowInd) > 0:
			for e in event.allBelowInd:
				if ((eventsLib.events[e].veamAges[-1] < AcceptableAge_MAX) and 
						(eventsLib.events[e].veamAges[-1] != -9999)):
					AcceptableAge_MAX = eventsLib.events[e].veamAges[-1]
		
		#Test for valid age range
		if AcceptableAge_MIN >= AcceptableAge_MAX:
			#bad_range = str(AcceptableAge_MAX-AcceptableAge_MIN)
			sys.stderr.write('\n	ERROR: NO ACCEPTABLE AGE RANGE FOR EVENT \'%s\'\n' % event.id)
			sys.stderr.write('AcceptableAge_MIN: %0.3f, AcceptableAge_MAX: %0.3\n' % (AcceptableAge_MIN, AcceptableAge_MAX))
			
		# This is where the fun begins
		mu    = event.ageModel[event.modelChoice]
		sigma = event.uncertModel[event.modelChoice]
		
		if sigma == 0: # If the date is historic or exact...
			event.veamAges[-1] = mu
			continue

		if sigma != 0:
			if mu == -9999:
				event.veamAges[-1] = np.random.uniform(AcceptableAge_MAX, AcceptableAge_MIN)
				continue
			else:
				# Convert minage and maxage to standard normal range because a, b are the standard deviations
				a = (AcceptableAge_MIN - mu) / sigma
				b = (AcceptableAge_MAX - mu) / sigma
				
				#Use Random Uniform if the area is very narrow. 
				#This is because of likelihood of bad result with truncnorm
				if (b-a) < 0.1:
					np.random.uniform(AcceptableAge_MAX, AcceptableAge_MIN)
				
				# Use truncated normal distribution to sample age, make sure it is greater than zero
				event.veamAges[-1] = truncnorm.rvs(a, b, loc=mu, scale=sigma)
				breakpt = 0 #break after 10
				if event.veamAges[-1] <=0:
					while ((event.veamAges[-1] <= 0) and (breakpt < 10)):
						event.veamAges[-1] = truncnorm.rvs(a, b, loc=mu, scale=sigma)
						breakpt += 1
				if np.isinf(event.veamAges[-1]):
					while ((np.isinf(event.veamAges[-1])) and (breakpt < 10)):
						event.veamAges[-1] = truncnorm.rvs(a, b, loc=mu, scale=sigma)
						breakpt += 1
				
				#If the sample is too far along the tail, just throw the age model out!
				# and choose from a random uniform.
				if breakpt == 10:
					event.veamAges[-1] = np.random.uniform(AcceptableAge_MAX, AcceptableAge_MIN)

				#Appen
				if event.veamAges[-1] < AcceptableAge_MIN:
					sys.stderr.write('This might be a minage issue!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
					return -1

				if event.veamAges[-1] > AcceptableAge_MAX:
					sys.stderr.write('This might be a maxage issue!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
					return -1

	return 0

def veam_main(inputs, VEAMWin):
	#Inputs needs to be an instance of a Variable Class
	
	#Define global minimum/maximum model age
	
	global GLOBAL_MINAGE
	global GLOBAL_MAXAGE
	
	start = time.time()
	
	sys.stdout.write('\nBeginning VEAM Simulation\n')
	
	GLOBAL_MINAGE = inputs.minAge
	GLOBAL_MAXAGE = inputs.maxAge
	#Define other variables
	ages_db_file = inputs.ageDB
	strat_db_file = inputs.stratDB
	use_mag = inputs.geoMag
	Rdt = inputs.res
	numruns = inputs.sims
	style = inputs.sorting
	
	
	errFlag, field = load_databases(strat_db_file,ages_db_file)
	if errFlag != 0:
		sys.stderr.write('Error while Loading Databases\n')
		return -1, None
	
	
	if style == 'ignore_strat':
		#If sorting is ignore strat, just remove the stratigraphy entries
		field.sortIgnoreStrat()
	
	for sim in range(numruns):
		if style == 'random':
			#Sort randomly
			field.sortRandom()
		elif style == 'most_contacts':
			#Sort by Number of Strat Contacts
			field.sortMostContacts()
		elif style == 'crater_age_uncertainty':
			#Sort by Age Uncertainty
			field.sortAgeUncertainty()
		elif style != 'ignore_strat':
			#No more sorting strategies left, error out.
			sys.stderr.write('Invalid sorting technique!')
			return -1, None
		
		#Get indices to all expanded strat relationships for fast min/max age finding later
		field.getStratIndices()
		#sys.stdout.write('\nEvents Sorted according to Age Uncertainty\n')
		
		ret = sample_ages(field)
		if ret != 0:
			sys.stdout.write('\nError in dating events!\n')
			return -1
		
		VEAMWin.veamProgress = int(((sim+1)*100.0)/numruns)
		
	
	elapsed = time.time() - start
	sys.stdout.write('     COMPLETED %d simulations in %0.4f seconds\n' % (numruns, elapsed))
	
	#Sort events alphabetically before printing out
	field.sortABC()
	return 0, field

### VEAM Variables and Structures ###

class Event():
	def __init__(self):
		self.id = ''
		self.ageModel    = [] # list of modeled ages
		self.uncertModel = [] # list of associated uncertainties
		self.veamAges    = [] # VEAM-simulated ages
		self.polarity    = None #Normal='N', Reverse='R'?
		self.stratAbove  = [] # ids of events that are immediately stratigraphically above
		self.stratBelow  = [] # ids of events that are immediately stratigraphically lower
		self.allAbove    = [] # ids of events indirectly above this event
		self.allBelow    = [] # ids of events indirectly below this event
		self.allAboveInd = [] # event library indices of allAbove
		self.allBelowInd = [] # event library indices of allBelow
		self.totalRels   = 0  # total number of stratigraphic relationships
		self.modelChoice = 0  # Which age model is being used in case of multiple models
		self.sortOrder   = 0  # Where is this event in line to be dated in VEAM?
		
	def addAgeModel(self, time, uncertainty):
		self.ageModel.append(float(time))
		self.uncertModel.append(float(uncertainty))
	
	def checkStratRels(self):
		#This is a naive check to see if there is a direct contridiction in the
		#Stratigraphy database.
		for rel in self.stratBelow:
			if rel in self.stratAbove:
				return False #if a strat relationship is in both upper and lower lists, that's an error
		return True #these lists are ok.
		
	def calcTotalRelCount(self):
		self.totalRels = len(self.stratAbove) + len(self.stratBelow)
	
	def chooseAgeModel(self):
		if len(self.ageModel) > 1:
			self.modelChoice = np.random.randint(len(self.ageModel))
	
class eventLibrary():
	def __init__(self):
		self.events = []
		self.tempEList = []

	def addEvent(self,name):
		### Create a new event
		e = Event()
		e.id = str(name)
		self.events.append(e)
		
	def checkAllStrat(self,relList):
		### Check Stratigraphy Web. Needs relationships, an Nx2 array of all 
		# Stratigraphic Relationships in StratDB file.
		#   If valid, return True, else return False
		for e in self.events:
			'''
			check = e.checkStratRels()
			if check == False:
				return False
			'''
			eID = e.id
			#print("Interrogating event:", eID)
			errorFlag = self.findProblemStrat(eID, relList, [eID]) 
			if len(errorFlag):
				sys.stderr.write("\n\nERROR: There is a stratigraphic loop with\n these vents:\n")
				for event in errorFlag:
					sys.stderr.write("  - "+event+"\n")
				return -1
		return 0
		
		#make a deeper check
		
		return True
	
	def findProblemStrat(self,curEvent,relList,eventsAbove):
		#This is similar to getEventsBelow, but specifically finds problematic 
		#Contridictions in the Stratigraphy Web.
		for r in relList:
			if curEvent == r[1]: #If the current event is in the upper position
				#print(r[0],r[1])
				#print ('event ', r[1], ' is above ', r[0])
				if r[0] in eventsAbove:
					#THIS MEANS THERE IS AN ERROR and it was JUST FIRST CAUGHT NOW
					#Crop the eventsabove list to exclude events outside the loop
					problemEvents = eventsAbove[eventsAbove.index(r[0]):]
					return problemEvents
				else:
					newEventsAbove = eventsAbove + [r[0]]
					result = self.findProblemStrat(r[0],relList,newEventsAbove)
					if len(result) != 0:
						return result
		return []
	
	def chooseAgeModels(self):
		### Choose random age models for events with multiple age models
		for e in self.events:
			e.chooseAgeModel()
			
	def getFullStratLists(self, maxStratRels):
		### Find all events above and below each event.
		# MaxStratRels should be the number of relationships in the strat DB
		for event in self.events:
			#Find events above current event
			self.tempEList = []
			ret = self.getEventsAbove(event.id, maxStratRels)
			if ret != 0:
				sys.stderr.write('Error in getEventsAbove\n')
				return -1
			event.allAbove = sorted(set(self.tempEList)) #Store unique values in allAbove
			#Find events below current event
			self.tempEList = []
			ret = self.getEventsBelow(event.id, maxStratRels)
			if ret != 0:
				sys.stderr.write('Error in getEventsBelow\n')
				return -1
			event.allBelow = sorted(set(self.tempEList)) #Store unique values in allBelow
		return 0
		
			
	def getEventsAbove(self, eID, maxAbove):
		### Recursively finds all events above an event and makes a temporary list of their ids
		# maxAbove should be the number of relationships in the strat DB
		# tempEList should start empty
		if len(self.tempEList) >= maxAbove:
			sys.stderr.write('Infinite Loop somewhere in strat relationships (going up)\n')
			return -1
		for event in self.events:
			if event.id == eID:   #Find event ID
				if len(event.stratAbove) > 0:
					for higher in event.stratAbove: #for all event ids in stratAbove
						self.tempEList.append(higher)  #append the eID to the temporary list
						ret = self.getEventsAbove(higher, maxAbove)    #RECURSIVE GET EVENTS ABOVE
						if ret == -1:
							return -1
				break
		return 0 #return when no longer recurring
			
	def getEventsBelow(self, eID, maxBelow):
		### Recursively finds all events below an event and makes a temporary list of their ids
		# maxAbove should be the number of relationships in the strat DB
		# tempEList should start empty
		if len(self.tempEList) >= maxBelow:
			sys.stderr.write('Infinite Loop somewhere in strat relationships (going down)\n')
			return -1
		for event in self.events:
			if event.id == eID:   #Find event ID
				if len(event.stratBelow) > 0:
					for lower in event.stratBelow: #for all event ids in stratAbove
						self.tempEList.append(lower)  #append the eID to the temporary list
						ret = self.getEventsBelow(lower, maxBelow)    #RECURSIVE GET EVENTS BELOW
						if ret == -1:
							return -1
				break
		return 0 #return when no longer recurring
	
	def sortABC(self):
		self.events.sort(key=operator.attrgetter('id'))
	
	def getStratIndices(self):
		### Get eventlibrary index for all upper and lower strat relationships
		for event in self.events:
			event.allAboveInd = []
			event.allBelowInd = []
			if len(event.allAbove) > 0:
				for rel in event.allAbove:
					for i,e in enumerate(self.events):
						if e.id == rel: 
							event.allAboveInd.append(i)
							break
			if len(event.allBelow) > 0:
				for rel in event.allBelow:
					for i,e in enumerate(self.events):
						if e.id == rel:
							event.allBelowInd.append(i)
							break
		return 0
		
	def sortRandom(self):
		### Sorts all events randomly
		np.random.shuffle(self.events)
		self.chooseAgeModels() # Select which age model to use for each event
	
	def sortMostContacts(self):
		### Sorts all events so events with more contacts are dated first.
		# Events with the same number of contacts get a random shuffle
		sortCount = 0 #Will be assigned in order 0 -> N for each event
		self.chooseAgeModels() # Select which age model to use for each event
		
		#Calculate number of total strat relationships for each event
		for e in self.events:
			e.calcTotalRelCount()
		
		#Make array of all total strat relationship counts
		contactCount = np.array([e.totalRels for e in self.events])
		#Go through all unique numbers of contacts with a sorted set
		for cs in sorted(set(contactCount)):
			#Find events that have the same number of strat contacts
			equalEvents = np.where(contactCount==cs)[0]
			#Give them a random shuffle
			np.random.shuffle(equalEvents)
			#Assign them a sort order
			for e in equalEvents:
				self.events[e].sortorder = sortCount
				sortCount += 1
		
		#Sort all events based on sort order
		self.events.sort(key=operator.attrgetter('sortorder'))
	
	def sortIgnoreStrat(self):
		### A bit of a misnomer, just gets rid of all stratigraphic relationships
		for e in self.events:
			e.stratAbove  = []
			e.stratBelow  = []
			
		self.chooseAgeModels() # Select which age model to use for each event
	
	def sortAgeUncertainty(self):
		### Sorts all events so events with least age model uncertainty are dated first.
		# Events with the same age uncertainty will get a random shuffle.
		# Events with more than one age model will have one chosen for them.
		sortCount = 0 #Will be assigned in order 0 -> N for each event
		self.chooseAgeModels() # Select which age model to use for each event
		
		for e in self.events:
			if e.uncertModel[0] < 0: #the uncertainty should only be -9999 if it's the only entry
				e.uncertModel[0] = UINT_MAX #change from -9999 to uintmax
		
		#Make an array of all the uncertainty values
		uncertainties = np.array([e.uncertModel[e.modelChoice] for e in self.events])
		for us in sorted(set(uncertainties)): #Go through all unique uncert vals
			equalEvents = np.where(uncertainties==us)[0] #All events with same uncert
			np.random.shuffle(equalEvents)               #Shuffle these events
			for e in equalEvents:                        #Assign each event a sort order
				self.events[e].sortorder = sortCount
				sortCount += 1
			if us == UINT_MAX:                           #Return unaged events to -9999
				for e in equalEvents:
					self.events[e].uncertModel[0] = -9999
		
		#Sort all events based on sort order
		self.events.sort(key=operator.attrgetter('sortorder'))
	
	
	
class variables():
	### Set of VEAM Variables	
	def __init__(self):
		
		self.minAge  = 0.0
		self.maxAge  = 0.0
		self.res     = 0.0
		self.ageUnit = ''
		self.ageDB   = ''
		self.stratDB = ''
		self.sims    = 0
		self.geoMag  = False
		self.igStrat = False
		self.sorting = ''
		
		#Update this if/when there are new or different styles of sorting by stratigraphy
		self.availSortStyles = ['random', 'crater_age_uncertainty', 'most_contacts', 'ignore_strat']
	
	def setMinAge(self,t):
		try:
			self.minAge = float(t)
		except ValueError:
			return False
		return True
		
	def setMaxAge(self,t):
		try:
			self.maxAge = float(t)
		except ValueError:
			return False
		return True
	
	def setAgeUnit(self,t):
		self.ageUnit = str(t)
		
	def setAgeDB(self,db):
		self.ageDB = str(db)

	def setStratDB(self,db):
		self.stratDB = str(db)
	
	def setResolution(self,t):
		try:
			self.res = float(t)
		except ValueError:
			return False
		return True
		
	def setSimulations(self,s):
		try:
			self.sims = int(s)
		except ValueError:
			return False
		return True
		
	def setGeoMag(self,b):
		if b=='True':
			self.geoMag = True
		else:
			self.geoMag  = False
		return True
		
	def setSorting(self,style):
		if style in self.availSortStyles:
			self.sorting = style
			return True
		else:
			return False
		
### GUI WINDOWS ###
class PlotCanvas(FigureCanvas):
	def __init__(self, parent=None, width=5, height=4, dpi=100, allCDFs=None,time=None, unit=''):
		fig = Figure(figsize=(width, height), dpi=dpi)
		
		FigureCanvas.__init__(self, fig)
		self.setParent(parent)
 
		FigureCanvas.setSizePolicy(self,
				QtWidgets.QSizePolicy.Expanding,
				QtWidgets.QSizePolicy.Expanding)
		FigureCanvas.updateGeometry(self)
		
		plotnum = min(len(allCDFs),1000)
		
		self.ax = self.figure.add_subplot(111)
		for i in range(plotnum):
			self.ax.plot(time,allCDFs[i],c='0.95')
		self.ax.plot(time,allCDFs[0],c='0.95',label=('%d Individual Simulations' % plotnum))
		self.ax.plot(time,np.mean(allCDFs, axis=0),c='r',label='Mean Cumulative Events')
		self.ax.plot(time,np.percentile(allCDFs, 10, axis=0),c='r',ls='--',label='10%')
		self.ax.plot(time,np.percentile(allCDFs, 90, axis=0),c='r',ls='--',label='90%')
		self.ax.set_title('Cumulative Events with time')
		
		self.ax.set_xlim([min(time),max(time)])
		self.ax.invert_xaxis()
		self.ax.legend(loc='upper left')
		self.ax.set_title('Cumulative Events with time')
		self.ax.set_ylabel('Cumulative Event Count')
		self.ax.set_xlabel('%s before present' % unit)
		
		self.draw()

class resultsWindow(QtWidgets.QWidget):
	def __init__(self, parent=None, fieldResults=None, sims=0, ageMin=0, ageMax=0, outputHeader='', ageUnit=''):
		super(resultsWindow, self).__init__()
		self.field = fieldResults
		self.sims = int(sims)
		self.ageMin = float(ageMin)
		self.ageMax = float(ageMax)
		self.outputHeader = str(outputHeader)
		self.ageUnit = ageUnit
		self.init_ui()
	
	def init_ui(self):
		
		self.allCDFs = self.genCDF()
		self.graph = PlotCanvas(self, width=6, height=4,allCDFs=self.allCDFs,time=self.time,unit=self.ageUnit)
		
		self.bSaveFig  = QtWidgets.QPushButton('Save Figure', self)
		self.bSaveData = QtWidgets.QPushButton('Save Data', self)
		self.lblFig = QtWidgets.QLabel('',self)
		self.lblData = QtWidgets.QLabel('',self)
		btn_box = QtWidgets.QHBoxLayout()
		btn_box.addWidget(self.bSaveFig)
		btn_box.addWidget(self.bSaveData)
		lbl_box = QtWidgets.QHBoxLayout()
		lbl_box.addWidget(self.lblFig)
		lbl_box.addStretch()
		lbl_box.addWidget(self.lblData)
		all_box = QtWidgets.QVBoxLayout()
		all_box.addWidget(self.graph)
		all_box.addLayout(btn_box)
		all_box.addLayout(lbl_box)
		
		self.setLayout(all_box)
		
		self.bSaveData.clicked.connect(self.saveData)
		self.bSaveFig.clicked.connect(self.saveFig)
		
		self.setWindowTitle('VEAM Results')
		self.setGeometry(250, 200, 640, 500)
		
		self.show()
	
	def saveData(self):
		### Save results to a csv file
		# Columns are individual events
		# rows are individual VEAM Simulations
		
		filename = QtWidgets.QFileDialog.getSaveFileName(self, 'Save Results')
		
		try:
			with open(filename[0], 'w') as outf:
				outf.write(self.outputHeader)
				#Write event names
				for i,e in enumerate(self.field.events):
					if i!=0:
						outf.write(',')
					outf.write('%s' % e.id)
				outf.write('\n')
				#Write simulation dates
				for s in range(int(self.sims)):
					for i,e in enumerate(self.field.events):
						if i!=0:
							outf.write(',')
						outf.write('%0.3f' % e.veamAges[s])
					outf.write('\n')
		except EnvironmentError:
			self.lblData.setText('')
			return -1
		
		self.lblData.setText('Data Saved!')
		return 0
	
		
	def saveFig(self):
		filename = QtWidgets.QFileDialog.getSaveFileName(self, 'Save Configuration File')
		pngT = ['png','unknown']
		jpgT = ['jpg', 'jpeg']
		tifT = ['tif','tiff']
		
		if len(filename[0])==0:
			self.lblFig.setText('Figure Not Saved')
			return -1
			
		figType = filename[0].split('.')[-1]
		
		if figType in jpgT:
			self.graph.print_jpg(filename[0])
		elif figType in pngT:
			self.graph.print_png(filename[0])
		elif figType in tifT:
			self.graph.print_tif(filename[0])
		else:
			self.graph.print_png(filename[0]+'.png')
				
		
		self.lblFig.setText('Figure Saved!')
	
	def cumFunc(self, data, time):
		### Create a CDF Function given a list of event times (data) in time bins (time)
		events = len(data)                #returns (length,width) of data
		simulation_bins = range(events+1) #histogram bins for the cumulative functions  
		CDF_func = np.zeros(len(time))    #Empty Mean Cumulative Function list
		
		for i,t in enumerate(time): #For each time bin
			#Time Histogram...
			histogram_t = np.histogram(np.where(data>t)[0],bins=simulation_bins)[0]
			#Average of events in MC solutions with time (CDF)
			CDF_func[i] = np.mean(histogram_t) * events
		
		return CDF_func
		
	def genCDF(self):
		### Generate CDFs for all simulations
		# This routine comes from cdf2rr.py
		self.ageStep = (self.ageMax - self.ageMin)/50.0
		self.time = np.arange(self.ageMin,self.ageMax,self.ageStep)
		
		sims = len(self.field.events[0].veamAges)
		
		### CDF for each Simulation
		cdfEach = np.zeros(( sims,len(self.time) ))
		simAges = np.zeros(len(self.field.events))
		
		for s in range(sims):
			#Get ages for this simulation
			for i,e in enumerate(self.field.events):
				simAges[i] = float(e.veamAges[s])
			#Get CDF for this simulation
			cdfEach[s] = self.cumFunc(simAges,self.time)
		
		return cdfEach



class progressWindow(QtWidgets.QWidget):
	def __init__(self, parent=None):
		super(progressWindow, self).__init__()
		self.init_ui()
	
	def init_ui(self):
		self.progress = QtWidgets.QProgressBar(self)
		
		self.btn = QtWidgets.QPushButton("Results", self)
		self.btn.setEnabled(False)
		
		h_box = QtWidgets.QHBoxLayout() #Vertical box
		h_box.addWidget(self.progress)
		h_box.addWidget(self.btn)
		
		v_box = QtWidgets.QVBoxLayout()
		v_box.addWidget(QtWidgets.QLabel('Running Veam...'))
		v_box.addLayout(h_box)
		
		self.setLayout(v_box)
		self.setGeometry(250, 200, 400, 10)
		
		self.btn.clicked.connect(self.getResults)
		
		self.setWindowTitle('Running VEAM...')
		self.show()
		
	def getResults(self):
			#Open Results Window, close this window
			
			self.close()
	

class mainWindow(QtWidgets.QWidget):
	
	def __init__(self, parent=None):
		super(mainWindow, self).__init__()
		self.init_ui()
		
	def init_ui(self):
		self.chkct = 0
		
		### The Banner Logo for VEAM
		banner = QtWidgets.QHBoxLayout()
		self.bannerLabel = QtWidgets.QLabel()
		self.bannerLabel.setPixmap(QtGui.QPixmap('greeley_cone_banner_600px.png').scaled(600,
																													600, QtCore.Qt.KeepAspectRatio))
		banner.addStretch()
		banner.addWidget(self.bannerLabel)
		banner.addStretch()
		
		self.lblGreeley = QtWidgets.QLabel('<i>Banner Photograph by Ron Greeley</i>')
		hCredit = QtWidgets.QHBoxLayout() #Young-Old Hor Box
		hCredit.addStretch()
		hCredit.addWidget(self.lblGreeley)
		
		
		### Age Ranges
		self.lblYoung = QtWidgets.QLabel('Minimum Possible Age:')
		self.lblOld   = QtWidgets.QLabel('Maximum Possible Age:')
		self.lblRes   = QtWidgets.QLabel('Temporal Resolution:')
		self.lblUnit  = QtWidgets.QLabel('Age Unit:')
		self.leYoung  = QtWidgets.QLineEdit()
		self.leOld    = QtWidgets.QLineEdit()
		self.leRes    = QtWidgets.QLineEdit()
		self.cbxAgeUnits = QtWidgets.QComboBox()
		#self.optYoung = QtWidgets.QComboBox
		
		self.cbxAgeUnits.addItems(('Ga', 'Ma', 'ka', 'a'))
		
		hYoung = QtWidgets.QHBoxLayout() #Young Hor Entry Box
		hYoung.addWidget(self.leYoung)
		hOld = QtWidgets.QHBoxLayout() #Old Hor Entry Box
		hOld.addWidget(self.leOld)
		gAges = QtWidgets.QGridLayout() #Young-Old Grid
		gAges.setVerticalSpacing(2)
		gAges.setHorizontalSpacing(100)
		gAges.addWidget(self.lblYoung, 1, 1)
		gAges.addWidget(self.lblOld, 1, 2)
		gAges.addLayout(hYoung, 2, 1)
		gAges.addLayout(hOld, 2, 2)
		gAges.addWidget(self.lblRes, 3, 1)
		gAges.addWidget(self.lblUnit, 3, 2)
		gAges.addWidget(self.leRes, 4, 1)
		gAges.addWidget(self.cbxAgeUnits, 4, 2)
		'''
		hRes = QtWidgets.QHBoxLayout() #Age Resolution Hor Box
		hRes.addWidget(self.lblRes)
		hRes.addWidget(self.leRes)
		hRes.addStretch()
		hRes.addWidget(self.lblUnit)
		hRes.addWidget(self.cbxRes)
		hRes.addStretch()'''
		
		vAges = QtWidgets.QVBoxLayout() #Total Age Range Box Layout
		vAges.addLayout(gAges)
		#vAges.addLayout(hRes)
		
		gpAges = QtWidgets.QGroupBox('Potential Volcanic Age Range') #Age Span groupbox
		gpAges.setLayout(vAges)
		
		### File Paths
		self.lblAgeDB   = QtWidgets.QLabel('Age Database Path ')
		self.lblStratDB = QtWidgets.QLabel('Stratigraphy Database Path ')
		self.leAgeDB    = QtWidgets.QLineEdit()
		self.leStratDB  = QtWidgets.QLineEdit()
		self.bAgeDB     = QtWidgets.QPushButton('...')
		self.bStratDB   = QtWidgets.QPushButton('...')
		
		gPath = QtWidgets.QGridLayout() #File Paths group Layout
		gPath.setSpacing(2)
		gPath.addWidget(self.lblAgeDB, 1, 0)
		gPath.addWidget(self.leAgeDB, 1, 1)
		gPath.addWidget(self.bAgeDB, 1, 2)
		gPath.addWidget(self.lblStratDB, 2, 0)
		gPath.addWidget(self.leStratDB, 2, 1)
		gPath.addWidget(self.bStratDB, 2, 2)
		
		
		gpDBs = QtWidgets.QGroupBox('Databases') #Databases groupbox
		gpDBs.setLayout(gPath)
		
		### Options
		self.lblSims     = QtWidgets.QLabel('Number of Simulations')
		self.leSims      = QtWidgets.QLineEdit()
		self.ckGeoMag    = QtWidgets.QCheckBox('Use Geomagnetic Data?')
		self.rdIgStrat   = QtWidgets.QRadioButton('Ignore Stratigraphy')
		self.rdSRandom   = QtWidgets.QRadioButton('Random')
		self.rdSContact  = QtWidgets.QRadioButton('Most Stratigraphic Contacts')
		self.rdSUncert   = QtWidgets.QRadioButton('Least Age Uncertainty')
		
		self.ckGeoMag.setDisabled(True)
		hSimsChks = QtWidgets.QHBoxLayout() #sims, geomag, strat use box
		hSimsChks.addWidget(self.lblSims)
		hSimsChks.addWidget(self.leSims)
		hSimsChks.addStretch()
		hSimsChks.addWidget(self.ckGeoMag)
		gSort =  QtWidgets.QGridLayout()
		gSort.addWidget(self.rdSRandom,1,1)
		gSort.addWidget(self.rdSContact,1,2)
		gSort.addWidget(self.rdSUncert,2,1)
		gSort.addWidget(self.rdIgStrat,2,2)
		vSorting = QtWidgets.QVBoxLayout()
		vSorting.addLayout(gSort)
		self.gpSort = QtWidgets.QGroupBox('Sorting Style') #Options groupbox
		self.gpSort.setStyleSheet('QGroupBox {font-weight: bold;}')
		self.gpSort.setLayout(vSorting)
		vOptions = QtWidgets.QVBoxLayout() #Options Box Layout
		vOptions.addLayout(hSimsChks)
		vOptions.addWidget(self.gpSort)
		
		gpOptions = QtWidgets.QGroupBox('Options') #Options groupbox
		gpOptions.setLayout(vOptions)
		
		### Submit Button
		self.bSubmit = QtWidgets.QPushButton('Run VEAM')
		self.bCheck  = QtWidgets.QPushButton('Check Inputs')
		self.bSaveCfg  = QtWidgets.QPushButton('Save Configuration')
		self.bLoadCfg  = QtWidgets.QPushButton('Load Configuration')
		self.lblSubmit = QtWidgets.QLabel()
		self.lblWarn   = QtWidgets.QLabel()
		
		
		self.bSubmit.setStyleSheet('QPushButton {background-color: skyblue; font-weight: bold;}')
		self.lblSubmit.setStyleSheet('QLabel {font-weight: bold; color: red; font-style: italic;}')
		self.lblWarn.setStyleSheet('QLabel {font-weight: bold; color: orange; font-style: italic;}')
		hSubmit = QtWidgets.QHBoxLayout() #Submit Buttons Layout
		hSubmit.addWidget(self.bSubmit)
		hSubmit.addStretch()
		hSubmit.addWidget(self.bCheck)
		hSubmit.addWidget(self.bSaveCfg)
		hSubmit.addWidget(self.bLoadCfg)
		hWarnings = QtWidgets.QHBoxLayout() #Warning Flags Layout
		hWarnings.addWidget(self.lblSubmit)
		hWarnings.addStretch()
		hWarnings.addWidget(self.lblWarn)
		
		### Assemble in vertical box
		v_box = QtWidgets.QVBoxLayout() #Vertical box for all main window
		v_box.addLayout(banner)
		v_box.addWidget(gpAges)
		v_box.addWidget(gpDBs)
		v_box.addWidget(gpOptions)
		v_box.addLayout(hSubmit)
		v_box.addLayout(hWarnings)
		v_box.addStretch()
		v_box.addLayout(hCredit)
		
		
		self.setLayout(v_box)
		self.setWindowTitle('VEAM | Volcanic Event Age Model')
		self.setGeometry(200, 100, 600, 400)
		
		### Actions
		
		self.bSaveCfg.clicked.connect(self.saveCfg)
		self.bLoadCfg.clicked.connect(self.loadCfg)
		self.bCheck.clicked.connect(lambda: self.checkCfg())
		self.bSubmit.clicked.connect(self.runVeam)
		
		self.bAgeDB.clicked.connect(lambda: self.open_path(self.leAgeDB))
		self.bStratDB.clicked.connect(lambda: self.open_path(self.leStratDB))
		
		self.show()
	
	def results(self):
		### Prepare a header for outputting ascii data and open the results window
		
			self.outputheader =  'Parameters\n'
			self.outputheader = self.outputheader + ('  Age Database:          %s\n' % self.vAgeDB)
			self.outputheader = self.outputheader + ('  Stratigraphy Database: %s\n' % self.vStratDB)
			self.outputheader = self.outputheader + ('  Min Age: %s\t Max Age: %s\t Units: %s\n' 
											% (self.vMinAge, self.vMaxAge, self.unitStrAgeDB))
			self.outputheader = self.outputheader + ('  VEAM Simulations:      %s\n' % self.vSims)
			self.outputheader = self.outputheader + ('  Event Sorting Style:   %s\n' % self.vSorting)
			self.outputheader = self.outputheader + ('  Use Geomagnetic Data?  %s\n' % self.vGeoMag)
			self.lblSubmit.setText('Viewing Results')
			
			self.resultsWin = resultsWindow(fieldResults=self.field, sims=self.vSims,
																	ageMax=self.vMaxAge, ageMin=self.vMinAge,
																	outputHeader=self.outputheader, 
																	ageUnit= self.unitStrAgeDB, parent=self)
			#self.resultsWin.bSaveData.clicked.connect(self.saveResults)
		
	def runVeam(self):
		### Run a check function that formats the form data first!
		
		check = 0
		check = self.checkCfg()
		
		if check == 0:
			veamVars = variables()
			veamVars.setMaxAge(self.vMaxAge)
			veamVars.setMinAge(self.vMinAge)
			veamVars.setResolution(self.vResol)
			veamVars.setAgeUnit(self.unitStrAgeDB)
			veamVars.setAgeDB(self.vAgeDB)
			veamVars.setStratDB(self.vStratDB)
			veamVars.setSimulations(self.vSims)
			veamVars.setGeoMag(self.vGeoMag)
			veamVars.setSorting(self.vSorting)
			
			self.lblSubmit.setMinimumWidth(200)
			self.lblSubmit.setText('Running VEAM...')
			self.lblSubmit.repaint()
			self.lblWarn.setText('')
			self.lblWarn.repaint()
			self.bSubmit.setEnabled(False)
			self.bSubmit.repaint()
			
			
			#Open the results window!
			self.progressWin = progressWindow(parent=self)
			self.veamProgress = 0
			#Create Thread with que as the return variable host for error handling
			#This will allow VEAM to be multithreaded if desired
			que = queue.Queue()
			thr = threading.Thread(target= lambda q, var, win: q.put(veam_main(var, win)), args=(que, veamVars, self))
			thr.start() # runs VEAM as Thread
			
			#During VEAM run, allow the windows to update and be used if needed
			while thr.is_alive():
				QtCore.QCoreApplication.processEvents()
				self.progressWin.progress.setValue(self.veamProgress)
				time.sleep(0.1)
			
			thr.join() # wait till VEAM is done
			errorFlag, self.field = que.get() #Get the return from VEAM
			
			#errorFlag, self.field = veam_main(veamVars, self) #Run 1 instance of VEAM!
			
			self.bSubmit.setEnabled(True) #Re-enable the Run VEAM Button
			if errorFlag==0:
				self.lblSubmit.setText('VEAM Completed!!')
				sys.stdout.write('\nVEAM Completed!!\n\n')
				self.progressWin.progress.setValue(100)
				self.progressWin.btn.setEnabled(True)
				self.progressWin.btn.clicked.connect(self.results)
				
			else:
				self.lblSubmit.setText('VEAM had an error! Check console for info')
				self.progressWin.close()
	
	def saveCfg(self):
		
		### Run a check function that formats the form data first!
		check = 0
		check = self.checkCfg()
		
		### Select a file and save
		filename = QtWidgets.QFileDialog.getSaveFileName(self, 'Save Configuration File')
		try:
			with open(filename[0], 'w') as cfg:
				cfg.write('This is a configuration file for the VEAM.py program\n')
				cfg.write('VEAM.py opens this file to set variables and input data file names.\n')
				cfg.write('\n## Inputs ##\n')
				cfg.write('Maximum Field Age: %s\n' % self.vMaxAge)
				cfg.write('Minimum Field Age: %s\n' %self.vMinAge)
				cfg.write('Age Database: %s\n' % self.vAgeDB)
				cfg.write('Stratigraphy Database: %s\n' % self.vStratDB)
				cfg.write('\n')
				cfg.write('Use geomagnetic information: %s\n' % self.vGeoMag)
				cfg.write('Timestep: %s\n' % self.vResol)
				cfg.write('Simulations: %s\n' % self.vSims)
				cfg.write('Sorting Style: %s\n' % self.vSorting)
				cfg.write('Time Units: %s\n' % self.unitStrAgeDB)
		except EnvironmentError: # parent of IOError, OSError *and* WindowsError where available
			self.lblWarn.setText('Config not saved')
			check = -1
		
		if check > 0: #Print Result in Orange Warning Label
			self.lblWarn.setText('Warning: Some errors remain, saved anyway')
		elif check == 0:
			self.lblWarn.setText('Configuration Saved!')
		
	def checkCfg(self):
		### Check that the form is filled correctly
        ### ALSO THIS IS WHERE VARIABLES FROM THE GUI ARE ASSIGNED TO CODE IF THEY ARE VALID ###
		
		errStyle = 'QLineEdit {background-color: pink}'
		valStyle = 'QLineEdit {background-color: white}'
		errMsg   = []
		
		#Find ComboBox Units (Ga, Ma, ka, or a)
		self.unitStrAgeDB = str(self.cbxAgeUnits.currentText())
		
		try:                                  #Young lineEdit
			minA = float(self.leYoung.text()) #Call variable to float
			self.vMinAge  = '%0.6f' % minA        #Set string variable
			if minA >= 0:                         #Has to be a valid age
				self.leYoung.setStyleSheet(valStyle) #set valid style for lineEdit
			else:
				errMsg.append('Err: Min. Age < 0!')  #Add Error Flag
				self.leYoung.setStyleSheet(errStyle) #set error style for lineEdit
		except ValueError:                      #If Min Age isn't a valid float
			errMsg.append('Err: Min. Age Not Number!') #Add Error Flag
			self.leYoung.setStyleSheet(errStyle)  #set error style for lineEdit
			self.vMinAge  = str(self.leYoung.text()) #Set string variable still
		
		try:                                  #Old lineEdit
			maxA = float(self.leOld.text())
			self.vMaxAge  = '%0.6f' % maxA
			self.leOld.setStyleSheet(valStyle)
			if len(errMsg)==0:
				if maxA <= minA:
					self.leOld.setStyleSheet(errStyle)
					errMsg.append('Err: Max. Age < Min. Age!')
		except ValueError:
			errMsg.append('Err: Max. Age Not Number!')
			self.leOld.setStyleSheet(errStyle)
			self.vMaxAge  = str(self.leOld.text())
		
		try:                                     #Temporal Resolution
			res = float(self.leRes.text())
			self.vResol   = '%0.6f' % res
			if res > 0:  #Resolution must be positive
				self.leRes.setStyleSheet(valStyle)
				# Resolution must be smaller than half the range, if the range is valid
				if len(errMsg)==0:
					if ((res > ((maxA-minA)/2.0)) and (maxA-minA > 0)):
						errMsg.append('Err: Timestep Too Large!')
						self.leRes.setStyleSheet(errStyle)
			else:
				errMsg.append('Err: Timestep not positive!')
				self.leRes.setStyleSheet(errStyle)
		except ValueError:
			errMsg.append('Err: Timestep Not Number!')
			self.leRes.setStyleSheet(errStyle)
			self.vResol   = str(self.leRes.text())
		
		self.vAgeDB   = str(self.leAgeDB.text())     #Age Database
		if os.path.isfile(self.vAgeDB):  #Check if it's a file
			self.leAgeDB.setStyleSheet(valStyle)
		else:
			self.leAgeDB.setStyleSheet(errStyle)
			errMsg.append('Err: Age DB doesn\'t exist!')
		
		self.vStratDB   = str(self.leStratDB.text())   #Strat Database
		if os.path.isfile(self.vStratDB): 
			self.leStratDB.setStyleSheet(valStyle)
		else:
			self.leStratDB.setStyleSheet(errStyle)
			errMsg.append('Err: Strat DB doesn\'t exist!')
		
		if self.vStratDB == self.vAgeDB:            #Strat and Age Database Compare
			self.leAgeDB.setStyleSheet(errStyle)
			self.leStratDB.setStyleSheet(errStyle)
			errMsg.append('Err: Strat and Age DBs are identical!')
		
		try:                                        #Simulation Count
			s = float(self.leSims.text())
			self.vSims = '%d' % s
			self.leSims.setStyleSheet(valStyle)
			if s <= 0:
				self.leSims.setStyleSheet(errStyle)
				errMsg.append('Err: < 1 Simulation!')
		except ValueError:
			self.leSims.setStyleSheet(errStyle)
			self.vSims = str(self.leSims.text())
			errMsg.append('Err: Sim Count Not Number!')
	
		self.vSorting = 'NONE'                       #Sorting Style
		                                             #Set Box to valid
		self.gpSort.setStyleSheet('QGroupBox {background-color: none; font-weight:bold;}')
		if self.rdSRandom.isChecked():
			self.vSorting = 'random'
		elif self.rdSContact.isChecked():
			self.vSorting = 'most_contacts'
		elif self.rdSUncert.isChecked():
			self.vSorting = 'crater_age_uncertainty'
		elif self.rdIgStrat.isChecked():
			self.vSorting = 'ignore_strat'
		else:                                        #Set Box to Error
			self.gpSort.setStyleSheet('QGroupBox {background-color: pink; font-weight:bold;}')
			errMsg.append('Err: Sorting Style not set!')
		
		#Set GeoMag Check
		self.vGeoMag  = '%r' % self.ckGeoMag.isChecked()
		
		#Print Errors to Error Message Area
		if len(errMsg) > 1:
			self.lblSubmit.setText('%s (+%d more)' % (errMsg[0], (len(errMsg)-1)))
			self.lblWarn.setText('')
		elif len(errMsg) == 1:
			self.lblSubmit.setText(errMsg[0])
			self.lblWarn.setText('')
		else: 
			self.lblSubmit.setText('')
			self.lblWarn.setText('No Errors Identified')
		
		return len(errMsg)
	
	def loadCfg(self):
		filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Database File Path')
		
		try: #Open the file and read the lines
			with open(filename[0], 'r') as cfg:
				rawcfg = cfg.readlines()
		except EnvironmentError: # parent of IOError, OSError *and* WindowsError where available
			self.lblWarn.setText('Config not loaded')
			return -1
		
		inputs = [line.rstrip() for line in rawcfg] # Strip end of line whitespace
		try:
			l = inputs.index('## Inputs ##') # Where does the line ## Inputs: ## occur
		except ValueError:
			self.lblWarn.setText('Could Not Load Configuration! (Missing inputs header)')
			return -1
		
		self.leYoung.setText(inputs[l+2].split(': ')[1])
		self.leOld.setText(inputs[l+1].split(': ')[1])
		self.leAgeDB.setText(inputs[l+3].split(': ')[1])
		self.leStratDB.setText(inputs[l+4].split(': ')[1])
		self.leRes.setText(inputs[l+7].split(': ')[1])
		self.leSims.setText(inputs[l+8].split(': ')[1])
		self.ckGeoMag.setChecked(inputs[l+6].split(': ')[1]=='True')
		style = str(inputs[l+9].split(': ')[1])
		self.rdSContact.setChecked(style=='most_contacts')
		self.rdSRandom.setChecked(style=='random')
		self.rdSUncert.setChecked(style=='crater_age_uncertainty')
		self.rdIgStrat.setChecked(style=='ignore_strat')
		
		#Set Unit Combo Boxes to the unit that everything is saved in
		try:
			self.unitStrAgeDB = str(inputs[l+10].split(': ')[1])
			self.cbxAgeUnits.setCurrentText(self.unitStrAgeDB)
		except IndexError: #If no unit is found, just assume 'Ma'
			self.cbxAgeUnits.setCurrentText('Ma')
		
		self.checkCfg()
		
		self.lblWarn.setText('Configuration Loaded!')
		
	def open_path(self,le):
		filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Database File Path')
		le.setText(filename[0])


def main():
	### Define QT App Instance
	app = QtCore.QCoreApplication.instance()
	if app is None:
		app = QtWidgets.QApplication(sys.argv)

	### Create a Window
	VEAM_window = mainWindow()
	
	#Exit upon closed window
	sys.exit(app.exec_())

if __name__ == '__main__':
	main()