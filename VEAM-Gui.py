#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 16:42:04 2017

@author: jaricha4
"""

import sys, os, time, operator
from PyQt5 import QtWidgets, QtCore, QtGui
import numpy as np
from scipy.stats import norm, truncnorm
from scipy.interpolate import interp1d
#import pylab as plt
import random

global start
global UINT_MAX
UINT_MAX = np.iinfo('uint32').max
	
def minage_finder_debug(eventID,relDB,ageDB,statement_min):
	#Find youngest possible age for event, based on stratigraphy
	#Look in relationship db for instances where eventID is in the OLDER column
	#If age of YOUNGER event is dated, (is not -9999)
	  #test to see if it is the OLDEST age seen so far
	#For each YOUNGER event, send through this function recursively
	#If no instances of eventID seen in OLDER column, return GLOBAL_MINAGE
	
	minage = GLOBAL_MINAGE
	#statement = []
	
	for relationship in relDB:
		if eventID==relationship[0]:
			#print 'eventID, relationship:', eventID, relationship[0]
			statement_min.append('%s:%s %f' % (eventID,relationship[1],ageDB[relationship[1]]))
			#print eventID+" is older than "+relationship[1]
			#if younger vent has been aged and is older than current minimum age:
			if ageDB[relationship[1]] != -9999 and ageDB[relationship[1]] > minage:
				#reassign minimum age to this vent's age.
				minage = ageDB[relationship[1]]
				
			#Recursively call minage_finder for each younger vent
			minage_indirect, statement_min = minage_finder_debug(relationship[1],relDB,ageDB,statement_min)
			#If indirectly younger vent is aged and is older than cur. min. age:
			if minage_indirect > minage:
				#reassign minimum age to indirect minimun age.
				minage = minage_indirect
	return minage, statement_min
	

def maxage_finder_debug(eventID,relDB,ageDB,statement_max):
	#Find oldest possible age for event, based on stratigraphy
	#Look in relationship db for instances where eventID is in the YOUNGER column
	#If age of OLDER event is dated, (is not -9999)
	  #test to see if it is the YOUNGEST age seen so far
	#For each OLDER event, send through this function recursively
	#If no instances of eventID seen in OLDER column, return GLOBAL_MAXAGE.
	
	maxage = GLOBAL_MAXAGE
	#statement = []
		
	for relationship in relDB:
		if eventID==relationship[1]:
			statement_max.append('%s:%s %f' % (eventID,relationship[0],ageDB[relationship[0]]))
			#print eventID+" is younger than "+relationship[0]
			#if older vent has been aged and is younger than current maximum age:
			if ageDB[relationship[0]] != -9999 and ageDB[relationship[0]] < maxage:
				#reassign minimum age to this vent's age.
				maxage = ageDB[relationship[0]]
				
			#Recursively call maxage_finder for each older vent
			maxage_indirect, statement_max =  maxage_finder_debug(relationship[0],relDB,ageDB,statement_max)
			#If indirectly older vent is aged and is younger than cur. max. age:
			if maxage_indirect < maxage:
				#reassign maximum age to indirect maximun age.
				maxage = maxage_indirect

	return maxage, statement_max

#Two functions, which are recursive, to find oldest and youngest bin_ding ages in a list of partially aged events
'''
def minage_finder(eventID,relDB,ageDB):
	#Find youngest possible age for event, based on stratigraphy
	#Look in relationship db for instances where eventID is in the OLDER column
	#If age of YOUNGER event is dated, (is not -9999)
	  #test to see if it is the OLDEST age seen so far
	#For each YOUNGER event, send through this function recursively
	#If no instances of eventID seen in OLDER column, return GLOBAL_MINAGE
	
	minage = GLOBAL_MINAGE
		
	for relationship in relDB:
		if eventID==relationship[0]:
			#print eventID+" is older than "+relationship[1]
			#if younger vent has been aged and is older than current minimum age:
			if ageDB[relationship[1]] != -9999 and ageDB[relationship[1]] > minage:
				#reassign minimum age to this vent's age.
				minage = ageDB[relationship[1]]
				
			#Recursively call minage_finder for each younger vent
			minage_indirect = minage_finder(relationship[1],relDB,ageDB)
			#If indirectly younger vent is aged and is older than cur. min. age:
			if minage_indirect > minage:
				#reassign minimum age to indirect minimun age.
				minage = minage_indirect

	return minage
	
def maxage_finder(eventID,relDB,ageDB):
	#Find oldest possible age for event, based on stratigraphy
	#Look in relationship db for instances where eventID is in the YOUNGER column
	#If age of OLDER event is dated, (is not -9999)
	  #test to see if it is the YOUNGEST age seen so far
	#For each OLDER event, send through this function recursively
	#If no instances of eventID seen in OLDER column, return GLOBAL_MAXAGE.
	
	maxage = GLOBAL_MAXAGE
		
	for relationship in relDB:
		if eventID==relationship[1]:
			#print eventID+" is younger than "+relationship[0]
			#if older vent has been aged and is younger than current maximum age:
			if ageDB[relationship[0]] != -9999 and ageDB[relationship[0]] < maxage:
				#reassign minimum age to this vent's age.
				maxage = ageDB[relationship[0]]
				
			#Recursively call maxage_finder for each older vent
			maxage_indirect = maxage_finder(relationship[0],relDB,ageDB)
			#If indirectly older vent is aged and is younger than cur. max. age:
			if maxage_indirect < maxage:
				#reassign maximum age to indirect maximun age.
				maxage = maxage_indirect
	return maxage
'''

def minunit_finder(eventID,relDB):
	#Find youngest possible age for event, based on stratigraphy
	#Look in relationship db for instances where eventID is in the OLDER column
	#If age of YOUNGER event is dated, (is not -9999)
	  #test to see if it is the OLDEST age seen so far
	#For each YOUNGER event, send through this function recursively
	#If no instances of eventID seen in OLDER column, return GLOBAL_MINAGE
	bin_= []
	for relationship in relDB:
		if eventID==relationship[0]:
			#print eventID+" is older than "+relationship[1]
			bin_.append(relationship[1])
			#print 'bin_', bin_

			#print 'Looking up a relation for:', relationship[1]
			#if relationship[1] in relationships[:,0]:
			minunit_indirect = minunit_finder(relationship[1],relDB)
			#print minunit_indirect

	return relationship[0], 'is older than', relationship[1], bin_

def maxunit_finder(eventID,relDB):
	#Find oldest possible age for event, based on stratigraphy
	#Look in relationship db for instances where eventID is in the YOUNGER column
	#If age of OLDER event is dated, (is not -9999)
	  #test to see if it is the YOUNGEST age seen so far
	#For each OLDER event, send through this function recursively
	#If no instances of eventID seen in OLDER column, return GLOBAL_MAXAGE.
	bin_ = []
	for relationship in relDB:
		if eventID==relationship[1]:
			#print eventID+" is younger than "+relationship[0]
			bin_.append(relationship[0])
			#print 'bin_', bin_

			#if relationship[0] in relationships[:,1]:
			maxunit_indirect = maxunit_finder(relationship[0], relDB)
				#bin_.append(maxunit_indirect[2])
				#print 'bin_', bin_

	return relationship[1], 'is younger than', relationship[0], bin_

def load_databases(strat_db_file, ages_db_file):
	#load age database using genfromtxt into 2xN string matrix
	relationships = np.genfromtxt(strat_db_file,skip_header=1,delimiter=',',dtype=None)
	#load ages database using genfromtxt
	try:
		AgeUncertainty = np.genfromtxt(ages_db_file, skip_header=1,delimiter=',', dtype=None)
	except:
		ValueError
		print '\n ERROR: Check for extra spaces in the names of events in %s and try again\n' % ages_db_file
		quit()        
	'''
	print '\nThis is the age database\n', AgeUncertainty
	'''
	MultipleUnitCheck = []
	Ages = {}
	Uncertainty = {}
	Polarity = {}
	
	for i in range(len(AgeUncertainty)): 
		MultipleUnitCheck.append(AgeUncertainty[i][0])

	MultipleUnitCheck = MultipleUnitCheck
	MUCset = [entry for entry in set(MultipleUnitCheck)]

	#print 'MultipleUnitCheck:', MultipleUnitCheck
	#print 'MUCset:', MUCset
	
	if len(MultipleUnitCheck) == len(MUCset): # If there is only one age per event...
		for i in range(len(AgeUncertainty)):
			Ages[AgeUncertainty[i][0]] = float(AgeUncertainty[i][1])
			Uncertainty[AgeUncertainty[i][0]] = float(AgeUncertainty[i][2])
			Polarity[AgeUncertainty[i][0]] = AgeUncertainty[i][4].lower()
	if len(MultipleUnitCheck) != len(MUCset):
		for i in range(len(AgeUncertainty)):
			Ages[AgeUncertainty[i][0]] = float(AgeUncertainty[i][1])
			Uncertainty[AgeUncertainty[i][0]] = float(AgeUncertainty[i][2])
			Polarity[AgeUncertainty[i][0]] = AgeUncertainty[i][4]
			#print AgeUncertainty[i][0], '=', AgeUncertainty[i][1]
		#print Ages
		CountAges = np.array([MultipleUnitCheck.count(entry) for entry in MUCset])
		MultipleAges = np.array(Ages.keys())[np.where(CountAges > 1)[0]]
		MultipleUnitCount = np.array(MultipleUnitCheck)
		for event in MultipleAges:
			eventages = AgeUncertainty[np.where(event == MultipleUnitCount)]
			eventageslist = []
			eventuncerlist = []
			eventpollist = []
			for age in eventages:				
				#print event, '=', age[1], age[2] #AgeUncertainty[np.where(event == MultipleUnitCount)]
				#print event, '=', age[1], age[2], age[4]
				eventageslist.append(float(age[1]))
				eventuncerlist.append(float(age[2]))
				eventpollist.append(age[4].lower())
			Ages[event] = eventageslist
			Uncertainty[event] = eventuncerlist
			Polarity[event] = eventpollist
	'''
	print '\nEvent Polarity'
	'''
	for key in Polarity:
		'''
		print key, Polarity[key]
		'''
		if 'normal' in Polarity[key]:
			if 'reversed' in Polarity[key]:
				print '\nPolarity reported for %s is both normal and reversed. Please update database file.' % key
				quit()
			Polarity[key] = 'n'

		if 'reversed' in Polarity[key]:
			if 'normal' in Polarity[key]:
				print '\nPolarity reported for %s is both normal and reversed. Please update database file.' % key
				quit()
			Polarity[key] = 'r'

		if len(Polarity[key]) != 1:
			Polarity[key] = 'NA'

	#print 'These should be the same number!!!', len(Polarity), len(set(Ages.keys()))

	Order_in = []
	events = []
	for idx in range(len(AgeUncertainty)):
		Order_in.append(AgeUncertainty[idx][3])
		events.append(AgeUncertainty[idx][0])

	Order_in_set = np.sort(np.array(list(set(Order_in))))
	events = np.array(events)
	Order_list = []
	#for order in Order_in_set:
		#indicies = np.where()
		#print order, list(events[Order_in == order]), events, events[Order_in == order]
		#Order_list.append(list(events[Order_in == order]))

	return relationships, Ages, Uncertainty, Polarity, Order_list#, trueages
'''
def load_mag_timescale():
  mag_timescale_file = 'Geomag_timescale_improved.csv'
  magts = np.genfromtxt(mag_timescale_file, delimiter=',', dtype=str)

  # Generate arrays of top and base age of chron with respect to polarities
  top_col = np.where(magts[0] == 'top (Ma)')
  base_col = np.where(magts[0] == 'base (Ma)')
  tops = np.ravel(np.transpose(magts[3:-1])[top_col]).astype(float) # Uses 3:-1 indicies, first 3 lines are headers and last one is incomplete
  bases = np.ravel(np.transpose(magts[3:-1])[base_col]).astype(float) # Uses 3:-1 indicies, first 3 lines are headers and last one is incomplete

  MINAGE_CUTOFF = 0
  if np.min(np.where(tops >= GLOBAL_MINAGE)) != 0:
    MINAGE_CUTOFF = np.min(np.where(tops >= GLOBAL_MINAGE)) - 1

  #print MINAGE_CUTOFF

  MAXAGE_CUTOFF = len(bases)
  if np.max(np.where(bases <= GLOBAL_MAXAGE)) != len(bases):
    MAXAGE_CUTOFF = np.max(np.where(bases <= GLOBAL_MAXAGE)) + 4

  #print MAXAGE_CUTOFF

  tops =  tops[MINAGE_CUTOFF:MAXAGE_CUTOFF]
  bases = bases[MINAGE_CUTOFF:MAXAGE_CUTOFF]
  #print tops, GLOBAL_MINAGE, bases, GLOBAL_MAXAGE

  # Generate an array of normal vs reversed chrons by looking at the last letter of the subchron column
  subchron_col = np.where(magts[0] == 'subchron')
  polarities = np.ravel(np.transpose(magts[3:-1])[subchron_col]) # Uses 3:-1 indicies, first lines are headers and last one is incomplete
  polarities = polarities[MINAGE_CUTOFF:MAXAGE_CUTOFF]
  polarities = np.array([p[-1:] for p in polarities]) # Polarity is the last char of each string in polarities
  pol_n = np.where(polarities == 'n')[0] # Indicies where polarity is normal
  pol_r = np.where(polarities == 'r')[0] # Indicies where polarity is reversed
  #print pol_n, pol_r

  #print tops
  #print tops[pol_n]
  #print bases[pol_n]

  return pol_n, pol_r, tops, bases
'''
'''
def make_topdown_chartlist():
  olderrelations = {}
  # Make a dictionary of all olderrelations for each event
  for i in range(len(Ages)):
    #print '########################### Now working on unit', Ages.keys()[i]
    #bin_ = []
    key, string, value, bin_ = minunit_finder(Ages.keys()[i], relationships)
    #print 'This is a key', Ages.keys()[i], string, value, set(bin_)
    setlist = [item for item in set(bin_)]
    olderrelations[Ages.keys()[i]] = setlist

  #print 'Out of the loop'

  # Convert the olderrelations dictionary to two lists
  eventkeys = olderrelations.keys()
  eventvalues = olderrelations.values()
  higher = [] # This is a list of events that are stratigraphically higher
  # This is a while loop that will start at the highest stratigraphic order and work down
  topdown = {}
  while len(eventkeys) > 0:
    #print len(eventkeys)
    #This is the initial case
    bin_ = []
    for k,v in zip(eventkeys, eventvalues): 
      if len(v) == 0:
        #print k, 'is one of the top strat units'
        bin_.append(k)
        eventkeys.remove(k)
        eventvalues.remove(v)
        higher.append(k)
    topdown[0] = bin_
    #This is for all subsequent cases
    for case in range(1, len(eventkeys)): # Won't actually go through the loop this many times!
      bin_ = []
      higherbin_ = []
      for k,v in zip(eventkeys, eventvalues):
        #print '############################ Working on event:', k
        nocount = 0
        yescount = 0
        for ik in olderrelations[k]:
          if ik in higher:
            #print 'yes', ik, 'is in', higher
            yescount += 1
          if ik not in higher:
            #print 'no', ik, 'is not in', higher
            nocount += 1
        if nocount == 0 and yescount != 0:
          #print k, 'is a second level strat unit'
          bin_.append(k)
          eventkeys.remove(k)
          eventvalues.remove(v)
          higherbin_.append(k)
      for usedevents in higherbin_:
        higher.append(usedevents)
      topdown[case] = bin_

  good = [item for item in topdown.values() if len(item) > 0]
  for ig in range(len(good)):
    #print 'Strat level', ig, good[ig]
    continue
  return olderrelations, [item for item in topdown.values() if len(item) > 0] # Only return the dictionary entries with data


def MakeStratGraph(runID, SampledAges, error_at, error_min, error_max, pastevents):
  olderrelations, events = make_topdown_chartlist()
  #print 'Now printing olderrelations and events'
  #print olderrelations, '\n\n',  events

  # This bit of code is specific to the Mars dataset!
  # Sort the events, by stratigraphic level. This works because of the way we originally hand drew the graph
  sorteventsint = []
  for tier in events:
    #print tier
    entriesint = []
    for entry in tier:
      #print entry[-2:]
      entriesint.append(entry[-2:])
    entriesint = np.array(entriesint, dtype=int)
    sortentries = np.sort(entriesint)
    sorteventsint.append(sortentries)

  #print sorteventsint, '\n'
  sortevents = []
  for tier in sorteventsint:
    #print tier
    entries = []
    for entry in tier:
      #print entry
      if entry < 10:
        entries.append('V0%i' % (entry))
      else:
        entries.append('V%i' % (entry))
    sortevents.append(entries)

  #print sortevents

  #print len(sortevents)

  widths = [len(iw) for iw in sortevents]
  #print max(widths)
  maxwidth = max(widths)

  coordinates = {}
  for y, level in enumerate(sortevents[::-1]):
    # Here need to sort the level based upon what is in next level up
    # This is done by olderrelations['each unit in level']
    # Group them
  #  for idx, units in enumerate(level):
  #    print units, olderrelations[units]
    #print '##########################################################'
    #print 'On this level', sortevents[::-1][y]
    try:
      #print 'Next level', sortevents[::-1][y+1]
      nextlevel = sortevents[::-1][y+1]
    except:
      IndexError
    order = []
    #for relations in olderrelations[unit]:
      #if relations not in nextlevel:
        #print relations, 'not in', nextlevel
      #if relations in nextlevel:
        #continue
    for nextunit in nextlevel:
      #print 'looking for which olderrelations is this:', nextunit
      #print 'unit, nextunit, olderrelations[unit]', unit, nextunit, olderrelations[unit]
      for idx, unit in enumerate(level):
        #print unit, olderrelations[unit]
        for relations in olderrelations[unit]:
          if relations == nextunit:
            #print 'yes', unit, 'should be appended to order'
            if unit not in order:
              order.append(unit)
              #print 'order', order
    for check in level:
      if check not in order:
        order.append(check)
    #print 'order', order
    #print 'level', level
    #for x, column in enumerate(order):
      #plt.text(1.5*x+maxwidth/len(order), 2*y, column)
      #coordinates[column] = (1.5*x+maxwidth/len(order)-0.1, 2*y+0.25)
    import pylab as plt
    pastevents = np.array(pastevents)
    for x, column in enumerate(level):
      if SampledAges[column] != -9999:
        plt.text(1.5*x+maxwidth/len(level), 2*y-1, '%i: %s\n%0.3f' % (np.where(column == pastevents)[0], column, SampledAges[column]))
      if column == error_at:
        plt.text(1.5*x+maxwidth/len(level), 2*y, '%s' % (column), color='red')
      #else:
        #plt.text(1.5*x+maxwidth/len(level), 2*y-1, '%s\nNS' % (column))
      coordinates[column] = (1.5*x+maxwidth/len(level)-0.1, 2*y+0.25)


  #print coordinates
  for ids in coordinates.values():
    plt.scatter(ids[0], ids[1], c='black')

  for older, youngers in zip(olderrelations.keys(), olderrelations.values()):
    if len(youngers) > 0:
      for younger in youngers:
        print older, younger
        xs = [coordinates[older][0], coordinates[younger][0]]
        ys = [coordinates[older][1], coordinates[younger][1]]
        plt.plot(xs, ys, c='blue', alpha=0.5)

  #plt.xlim(-1,len(events))
  plt.ylim(-1,2*len(events))
  #plt.show()
  plt.title('Error at: %s\nBetween: %s and %s' % (error_at, error_min, error_max))
  print 'This is runID:', runID
  plt.savefig('Stratgraph_%i.png' % (runID))
  plt.close()
'''

'''
Implement ranking strategies and results hypotheses

1. Random

2. Top down, Highest stratigraphic order first
   Should push ages older
plt.show()
3. Bottom up, lowest first
   Push ages younger

4. Most neighbors to least neighbors
   Start in middle, do top, do bottom

5. Most neighbors to least
   Randomly choose from list where most neighbors has more instances and more likely to be selected

6. Crater age uncertainty
   Do the lowest first

7. Gaussian from head and tail of group, uniform random for all else
   From Chuck

8. No strat, just sample ages
   From Chuck

alt. Split V00,V09 from V10,V17
'''
# These functions generate a list of events for sample_ages() to iterate over
def random_events(Ages):
  '''
  This function generates a random list of events for sample_ages() to iterate over
  '''
  events = [item for item in Ages]
  np.random.shuffle(events)
  return events
'''
def make_topdown_list():
  olderrelations = {}
  # Make a dictionary of all olderrelations for each event
  for i in range(len(Ages)):
    #print '########################### Now working on unit', Ages.keys()[i]
    #bin_ = []
    key, string, value, bin_ = minunit_finder(Ages.keys()[i], relationships)
    #print 'This is a key', Ages.keys()[i], string, value, set(bin_)
    setlist = [item for item in set(bin_)]
    olderrelations[Ages.keys()[i]] = setlist

  # Convert the olderrelations dictionary to two lists
  eventkeys = olderrelations.keys()
  eventvalues = olderrelations.values()
  higher = [] # This is a list of events that are stratigraphically higher
  # This is a while loop that will start at the highest stratigraphic order and work down
  topdown = {}
  while len(eventkeys) > 0:
    #This is the initial case
    bin_ = []
    for k,v in zip(eventkeys, eventvalues): 
      if len(v) == 0:
        #print k, 'is one of the top strat units'
        bin_.append(k)
        eventkeys.remove(k)
        eventvalues.remove(v)
        higher.append(k)
    topdown[0] = bin_
    #This is for all subsequent cases
    for case in range(1, len(eventkeys)+1): # Won't actually go through the loop this many times!
      bin_ = []
      higherbin_ = []
      for k,v in zip(eventkeys, eventvalues):
        #print '############################ Working on event:', k
        #print 'higher', higher
        nocount = 0
        yescount = 0
        for i in olderrelations[k]:
          if i in higher:
            #print 'yes', i, 'is in', higher
            yescount += 1
          if i not in higher:
            #print 'no', i, 'is not in', higher
            nocount += 1
        if nocount == 0 and yescount != 0:
          #print k, 'is a lower level strat unit'
          bin_.append(k)
          eventkeys.remove(k)
          eventvalues.remove(v)
          higherbin_.append(k)
      for usedevents in higherbin_:
        higher.append(usedevents)
      topdown[case] = bin_
      #print eventkeys

  #good = [item for item in topdown.values() if len(item) > 0]
  #for i in range(len(good)):
    #print 'Strat level', i, good[i]

  #print [item for item in topdown.values() if len(item) > 0]
  return [item for item in topdown.values() if len(item) > 0] # Only return the dictionary entries with data
'''
'''
def topdown_events(a):
  #''
  #This function takes a list of lists, where the first list is the highest stratigraphic order and the last list is the lowest stratigraphic order.
  #''
  events = []
  for i in a:
    np.random.shuffle(i)
    for event in i:
      events.append(event)

  #print events
  return events

def make_bottomup_list():
  youngerrelations = {}
  for i in range(len(Ages)):
    #print 'Now working on unit', Ages.keys()[i]
    key, string, value, bin_ = maxunit_finder(Ages.keys()[i], relationships)
    #print 'This is a key', Ages.keys()[i], string, value, set(bin_)
    setlist = [item for item in set(bin_)]
    youngerrelations[Ages.keys()[i]] = setlist

  eventkeys = youngerrelations.keys()
  eventvalues = youngerrelations.values()
  higher = [] # This is a list of events that are stratigraphically higher

  # This is a while loop that will start at the highest stratigraphic order and work down
  bottomup = {}
  while len(eventkeys) > 0:
    #This is the initial case
    bin_ = []
    for k,v in zip(eventkeys, eventvalues): 
      if len(v) == 0:
        #print k, 'is one of the top strat units'
        bin_.append(k)
        eventkeys.remove(k)
        eventvalues.remove(v)
        higher.append(k)
    bottomup[0] = bin_
    #This is for all subsequent cases
    for case in range(1, len(eventkeys)+1): # Won't actually go through the loop this many times!
      bin_ = []
      higherbin_ = []
      for k,v in zip(eventkeys, eventvalues):
        #print '############################ Working on event:', k
        nocount = 0
        yescount = 0
        for i in youngerrelations[k]:
          if i in higher:
            #print 'yes', i, 'is in', higher
            yescount += 1
          if i not in higher:
            #print 'no', i, 'is not in', higher
            nocount += 1
        if nocount == 0 and yescount != 0:
          #print k, 'is a second level strat unit'
          bin_.append(k)
          eventkeys.remove(k)
          eventvalues.remove(v)
          higherbin_.append(k)
      for usedevents in higherbin_:
        higher.append(usedevents)
      bottomup[case] = bin_

  #good = [item for item in bottomup.values() if len(item) > 0]
  #for i in range(len(good)):
    #print 'Strat level', i, good[i]

  return [item for item in bottomup.values() if len(item) > 0]
'''
'''
def bottomup_events(a):
  #''
  #This function takes a list of lists, where the first list is the lowest stratigraphic order and the last list is the highest stratigraphic order.
  #''
  events = []
  for i in a:
    np.random.shuffle(i)
    for event in i:
      events.append(event)
  return events
'''

def most_contacts_events(Ages, relationships):
  '''
  This function looks at the Stratigraphic relationships database and counts how many instances there are of each unit in both columns. Every instance represents a neighbor. 
  The idea here is to sample_ages of units with the most neighbors first and then sample_ages of units with progressively fewer neighbors. 
  This scheme allows for the first unit's to be sampled according to their probability distribution function, 
  with little to no influence from the minage/maxage recursive scheme which may trim the probability distribution 
  function of subsequently sampled events.
  '''
  #print relationships[:,0]
  #print relationships[:,1]
  neighbors = {}
  for i in Ages:
    neighbors[i] = len(np.nonzero(relationships[:,0] == i)[0]) + len(np.nonzero(relationships[:,1] == i)[0])
  neighborcounts = neighbors.values()
  neighborcounts = np.array(neighborcounts)
  neighborkeys = neighbors.keys()
  neighborkeys = np.array(neighborkeys)
  #print neighborcounts[::-1]
  events = []
  for i in range(max(neighborcounts), 0, -1):
    neighborset = np.where(neighborcounts == i)[0]
    #print neighborsetrelationships
    #print neighborkeys[neighborset]
    new = neighborkeys[neighborset]
    np.random.shuffle(new)
    for j in new:
      events.append(j)
  agekeys = np.array(Ages.keys())
  for i in agekeys:
    if i not in events:
      events.append(i)
      #print i, 'not in events'
  #print 'events', events
  #print 'len Events and len Ages', len(events), len(Ages)
  #print events, neighborkeys, neighborcounts
  return events

'''
def outside_in_events(top1, bottom1):
  #''
 #This function generates two outputs, outer and inner. It generates two outputs because it feeds a list of events to two different age sampling algorithms. The outer is comprised of units at the top and bottom of the stratigraphic sequence, thus there are components top1 and bottom1 from topdown and bottomdown, respectively. The inner is comprised of units that are in the middle of the stratigraphic sequence. The outer and inner events are randomly shuffled so that each time the function is called, there is a different sequence of events to iterate over. 
  #''
  outer = top1 + bottom1
  inner = [item for item in Ages if item not in outer]
  np.random.shuffle(outer)
  np.random.shuffle(inner)
  return outer + inner
'''

'''
def most_contacts_list_events():
  #''
  #This algorithm is very similar to most_contacts_events(), but there is randomness. This randomness results from a heirarchial scheme which puts more instances of a unit based upon number of neighbors. The number of instances is governed by taking the integer of e^(number_of_neighbors - 1), such that units with only one neighbor will only have one instance in the list. The number of instances scales exponentially as units have more neighbors. This final list is then shuffled so that all units in the list are in random locations. A unit is chosen at random from the list and then removed. The purpose is that units with the most neighbors will most likely be sampled first, but not all of the time.
  #''
  #print relationships[:,0]
  #print relationships[:,1]
  neighbors = {}
  for i in Ages:
    # This one liner alleviates the need for loops, providing optimization
    neighbors[i] = len(np.nonzero(relationships[:,0] == i)[0]) + len(np.nonzero(relationships[:,1] == i)[0])
  neighborcounts = neighbors.values()         # Make a list of values from the dictionary
  neighborcounts = np.array(neighborcounts)   # Convert the list to np.array so that we can use np.where()
  neighborkeys = neighbors.keys()             # Make a list of keys from the dictionary
  neighborkeys = np.array(neighborkeys)       # Convert the list to np.array so that we can use np.where()
  manyevents = []
  for i in range(max(neighborcounts), 0, -1):
    neighborset = np.where(neighborcounts == i)[0]
    #print neighborset
    #print neighborkeys[neighborset]
    new = neighborkeys[neighborset]
    newer = list(new)*int(np.e**(i-1))
    for j in newer:
      manyevents.append(j)
  for i in Ages:
    if i not in manyevents:
      manyevents.append(i)
  np.random.shuffle(manyevents)
  #print 'set Manyevents same as Ages', len(set(manyevents)) == len(Ages)
  events = []
  #print manyevents, '\n'
  for i in range(len(Ages)-1):
    if len(manyevents) > 1:
      choice = random.choice(manyevents)
      #print choice
      manyevents = [item for item in manyevents if item != choice]
    events.append(choice)
  events.append(manyevents[0])
    #print choice, manyevents, '\n'
  #print 'events', events
  #print 'length of events and ages', len(events), len(Ages)
  return events
'''

def crater_age_uncertainty_events(Uncertainty):
  '''
  This function samples units with the lowest age uncertainty first. When multiple units have the same uncertainty, the units are randomly shuffled. Thus, each time the function is called, there is a different list of events to iterate over.
  '''
  uncertainkeys = np.array(Uncertainty.keys())
  uncertainvalues = np.array(Uncertainty.values(),dtype=object)
  uncertainset = [item for item in list(Uncertainty.values())]
  #print sorted(uncertainset)
  events = []
  for i in sorted(uncertainset):
    uncertainsetbit = np.where(uncertainvalues == i)[0]
    np.random.shuffle(uncertainsetbit)
    #print uncertainsetbit
    for j in uncertainsetbit:
      #print uncertainkeys[j]
      events.append(uncertainkeys[j])
  
  return events
'''
def user_defined_events():
  events = []
  for sublist in Order_list:
    random.shuffle(sublist)
    events.extend(sublist)
  
  return events
'''
'''
def key_stratigraphic_unit_events():
  #''
  #This function looks at the Stratigraphic relationships database and counts how many instances there are of each unit in both columns. Every instance represents a neighbor. The idea here is to sample_ages of units with the most neighbors first and then sample_ages of units with progressively fewer neighbors. This scheme allows for the first unit's to be sampled according to their probability distribution function, with little to no influence from the minage/maxage recursive scheme which may trim the probability distribution function
  #''
  print relationships[:,0]
  print relationships[:,1]
  neighbors = {}
  event_ageKey = {}
  # Count the number of neighbors
  for i in Ages:
    # Need if, elif and else depending on number of events in strat database due to algorithm
    if len(relationships) > 2:
      neighbors[i] = len(np.nonzero(relationships[:,0] == i)[0]) + len(np.nonzero(relationships[:,1] == i)[0])
    elif len(relationships) == 2:
      neighbors[i] = len(np.nonzero(relationships[0] == i)[0]) + len(np.nonzero(relationships[1] == i)[0])
    else: # no need to go further because no strat order, just shuffle all of the events and go on
      neighbors[i] = 0

  # Use hierarchy of most neighbors to fewest
  neighborcounts = np.array(neighbors.values()) # neighbor counts
  neighborkeys = np.array(neighbors.keys())     # name of event corresponding to count
  #print neighborcounts[::-1]
  events = [] # empty list to append list of events
  for i in range(max(neighborcounts), 0, -1):
    neighborset = np.where(neighborcounts == i)[0]
    #print '\n', i, 'neighbors'                                             # PRINT FOR SUMMARY
    #print 'VentID \t Uncertainty'                                          # PRINT FOR SUMMARY
    #print neighborset
    #print neighborkeys[neighborset], 'new'
    new = neighborkeys[neighborset]
    uncers = np.zeros(len(new))
    for ventidx, ventid in enumerate(new):
      #print ventid, '\t', Uncertainty[ventid]                              # PRINT FOR SUMMARY
      try:
        uncers[ventidx] = Uncertainty[ventid]
      except:
        ValueError
        ageKey = np.random.randint(0,len(Uncertainty[ventid])) # choose a random index
        event_ageKey[ventid] = ageKey # assign the index, ageKey, to a dict for the sample_ages function to know which to use based upon order
        #print ventid, ageKey
        uncers[ventidx] = Uncertainty[ventid][ageKey] # Use this uncertainty for sorting the events
    #print 'Unique set of uncertainty values for this group of events'      # PRINT FOR SUMMARY
    #print np.unique(uncers)                                                # PRINT FOR SUMMARY
    #print 'Events sorted by uncertainty'                                   # PRINT FOR SUMMARY
    for uncer in np.sort(np.unique(uncers)):
      #print np.where(uncers == uncer)[0]
      uncers_by_group = np.where(uncers == uncer)[0]
      np.random.shuffle(uncers_by_group)
      #print uncers_by_group
      for j in uncers_by_group:
        events.append(new[j])
        #print new[j], uncers[j]                                            # PRINT FOR SUMMARY

  # This section is necessary in case an event has multiple ages. Since events are sorted by uncertainty, need to account for which age/uncertainty pair is to be used for this model run. For example, let's say an event has a precise Ar/Ar date and relatively imprecise K/Ar date whose uncertainty is greater than some other event, then the order of this event with respect to another depends upon which date we use. If the code randomly chooses the Ar/Ar date to sample first, then we need to be sure to use that date in the age assignment algorithm because the K/Ar date will invalidate the Key_Stratigraphic_Unit event order.
  agekeys = np.array(Ages.keys())
  agevalues = np.zeros(len(agekeys))
  for idx, age in enumerate(Ages.values()): # iterate through the elements in Ages.values to look for multiple ages per event
    try:
      type(age) == float
      agevalues[idx] = age
      event_ageKey[Ages.keys()[idx]] = 0
#      if age == -9999:
#        agevalues[idx] = -9999
#        event_ageKey[Ages.keys()[idx]] = 0
#      else:
#        agevalues[idx] = 1
#        event_ageKey[Ages.keys()[idx]] = 0
    except:
      TypeError
      count = 0 # set a counter
      for age_n in age: # loop through the ages in the list
        if age_n == -9999: # if the age is -9999, add one to the counter
          count += 1
      if count == len(age): # if all of the ages in list are -9999, then agevalues[idx] = -9999 to get it to go to end for uniform random
        agevalues[idx] = -9999
        event_ageKey[Ages.keys()[idx]] = 0
      else:
        agevalues[idx] = 1

  for i in agekeys:
    if i not in events:
      print i, 'was not incorporated into event order list. Appending to events' # These events don't have strat relations
      events.append(i)
      #print i, 'not in events'

  if -9999 in Ages.values():
    #print 'agekeys, agevalues', agekeys, agevalues
    Move2End = agekeys[np.where(agevalues == -9999)[0]]
    np.random.shuffle(Move2End)
    for unknownage in Move2End:
      events.remove(unknownage)
      events.append(unknownage)
  #print events
  #print 'events', events
  #print 'len Events and len Ages', len(events), len(Ages)
  #print events, neighborkeys, neighborcounts
  #print events
  #quit()
  #events = ['V13', 'V25', 'V02', 'V15', 'V03', 'V21', 'V05', 'V23', 'V11', 'V20', 'V10', 'V07', 'V12', 'V04', 'V14', 'V06', 'V19', 'V16', 'V17', 'V28', 'V27', 'V00', 'V24', 'V22', 'V26', 'V01', 'V18', 'V08', 'V09'] # THIS KEEPS SAME ORDER FOR 2sigfig vs 3sigfig for Mars
  #print events
  #quit()
  #print event_ageKey
  return events, event_ageKey
'''

def sample_ages(events, relationships, runID, Ages, Uncertainty, event_ageKey=None):
  '''
  This function generates a dictionary of sampled ages. 
  It iterates over a list of events. 
  It finds an acceptable age range based upon stratigraphic relationships and SampledAges. 
  This acceptable age range then truncates the probability distribution function. 
  A cumulative distribution function is created so that a random number from 0-1 
    can be selected and the associated age sampled, from the cumulative distribution function. 
  In cases where minage/maxage from stratigraphy truncate the probability distribution function 
    to one of the tails where the probability of an event is zero, the code then chooses a 
    random uniform age between the two bounding units. 
  In this way, stratigraphy can override the radiometric date.
  '''
  SampledAges = {} # Need a for loop to create a new dictionary for SampledAges
  #print 'relationships', relationships
  for i in Ages:
    SampledAges[i] = -9999 # Initially set SampledAges to -9999
  pastevents = []
  statement_min = []
  statement_max = []
  
  for currentevent in events:
    #print 'sampling event %s' % currentevent
    # This section was used for debugging.
    #sys.stdout.write("\nevent %s\t" % currentevent)
    statement_min.append('##%s##' % currentevent)
    statement_max.append('##%s##' % currentevent)
    pastevents.append(currentevent)
    #####FIND EVENT AGE RANGE######
    AcceptableAge_MIN, statement_min = minage_finder_debug(currentevent,relationships,SampledAges,statement_min)
    #sys.stdout.write('i')
    #print 'passed minage--debug!!! (Minumum age: %0.2f)' % AcceptableAge_MIN
    AcceptableAge_MAX, statement_max = maxage_finder_debug(currentevent,relationships,SampledAges,statement_max)
    #sys.stdout.write('a')
    #print 'passed maxage--debug!!! (Maximum age: %0.2f)' % AcceptableAge_MAX
    #AcceptableAge_MIN = minage_finder(currentevent,relationships,SampledAges)
    #AcceptableAge_MAX = maxage_finder(currentevent,relationships,SampledAges)
    #Test for valid age range
    if AcceptableAge_MIN >= AcceptableAge_MAX:
      bad_range = str(AcceptableAge_MAX-AcceptableAge_MIN)
      print "\n  ERROR: NO ACCEPTABLE AGE RANGE FOR EVENT ("+bad_range+" yrs)"+currentevent
      print AcceptableAge_MIN, AcceptableAge_MAX
      SampledAgeskeys = np.array(SampledAges.keys())
      SampledAgesvalues = np.array(SampledAges.values())
      print "minage at unit", SampledAgeskeys[np.where(SampledAgesvalues == AcceptableAge_MIN)], "age", SampledAgesvalues[np.where(SampledAgesvalues == AcceptableAge_MIN)]
      print "maxage at unit", SampledAgeskeys[np.where(SampledAgesvalues == AcceptableAge_MAX)], "age", SampledAgesvalues[np.where(SampledAgesvalues == AcceptableAge_MAX)]
      print 'This is the runID in sampleages():', runID
      print pastevents
      print '\nUnit is older than Unit Age'
      for entry in statement_min:
        print entry
      print pastevents
      print '\nUnit is younger than Unit Age'
      for entry in statement_max:
        print entry
      #print 'Statement_min', statement_min #debug
      #print 'Statement_max', statement_max #debug
      #print 'pastevents', pastevents #debug
      #print 'SampledAges', SampledAges
      error_at = currentevent
      error_min = SampledAgeskeys[np.where(SampledAgesvalues == AcceptableAge_MIN)]
      error_max = SampledAgeskeys[np.where(SampledAgesvalues == AcceptableAge_MAX)]
      #MakeStratGraph(runID, SampledAges, error_at, error_min, error_max, pastevents)
      break
    rangestr = "%0.2f-%0.2f" % (AcceptableAge_MIN,AcceptableAge_MAX)

    # This is where the fun begins
    if event_ageKey == None: # If we didn't use key_stratigraphic_unit for event sorting
      #sys.stdout.write('1')
      try:
        len(Ages[currentevent]) > 1 # If there is more than one age/uncertainty reported per event, randomly choose one.
        agechoice = np.random.randint(0, len(Ages[currentevent]))
        mu = Ages[currentevent][agechoice]
        sigma = Uncertainty[currentevent][agechoice]
        #sys.stdout.write('a')
      except TypeError:
        mu = Ages[currentevent]
        sigma = Uncertainty[currentevent]
        #sys.stdout.write('b')

    if event_ageKey != None:
      #sys.stdout.write('2')
      ageKey = event_ageKey[currentevent]
      #print ageKey
      if ageKey == 0: # if there is only one age
        try:
          len(Ages[currentevent]) > 1
          mu = Ages[currentevent][ageKey]
          sigma = Uncertainty[currentevent][ageKey]
          #print 'mu, sigma', mu, sigma
        except:
          TypeError
          mu = Ages[currentevent]
          sigma = Uncertainty[currentevent]
          #print 'mu, sigma', mu, sigma
      if ageKey > 0:
        mu = Ages[currentevent][ageKey]
        sigma = Uncertainty[currentevent][ageKey]

    if sigma == 0: # If the date is historic or exact...
      #sys.stdout.write('3')
      SampledAges[currentevent] = mu
      continue

    # Use mag age assignment routine
    # Set acceptable min/max age as the stratigraphic min/max age
    '''
    if use_mag == True and sigma == -9999: # If there is paleomag but not radiometric age
      #print currentevent, Polarity[currentevent]
      minage_strat = AcceptableAge_MIN
      maxage_strat = AcceptableAge_MAX
      #print 'Minage:', minage_strat, 'Maxage:', maxage_strat

      if Polarity[currentevent] == 'NA': # If paleomag, no radiometric age, polarity is 'NA'
        minage = AcceptableAge_MIN
        maxage = AcceptableAge_MAX
        SampledAges[currentevent] = np.random.uniform(maxage, minage)
        continue
        #print 'min, sampled, max', minage,SampledAges[currentevent], maxage

      range_strat = np.arange(minage_strat, maxage_strat, Rdt) # Timeline
      if Polarity[currentevent] == 'n': # If paleomag, no radiometric age, polarity is 'n'
        mag_function = mag_interpolate_n(range_strat)
        #plt.plot(range_strat, mag_function) #Polarities, 0 false, 1 true

        cdf = np.cumsum(mag_function)

        #''
        #        print cdf[-1]
        #        if cdf[-1] == 0:
        #          plt.title('%s %s\nmu: %0.2f sd: %0.2f' % (currentevent, Polarity[currentevent], mu, sigma))
        #          plt.ylim(-1, 2)
        #          plt.show()
        #          plt.close()
        #'

        if cdf[-1] == 0:
          print '\n', currentevent
          print 'range_strat', range_strat
          print 'mag_function', mag_function
          print 'pdf', pdf

        cdf /= np.max(cdf) # Normalize to one

        #plt.plot(range_strat, cdf)
        #plt.title('From mag only, no radiometric age')
        #plt.show()
        #plt.close()
          
        #print np.min(np.where( cdf >= np.random.uniform() )) * dt + minage
        SampledAges[currentevent] = np.min(np.where( cdf >= np.random.uniform() )) * Rdt + minage_strat
        #plt.axvline(SampledAges[currentevent], color='black', lw=2)

        #plt.title('%s %s\nmu: %0.2f sd: %0.2f' % (currentevent, Polarity[currentevent], mu, sigma))
        #plt.ylim(-1, 2)
        #plt.show()
        #plt.close()

      if Polarity[currentevent] == 'r':
        mag_function = mag_interpolate_r(range_strat)
        #plt.plot(range_strat, mag_function, lw=2, label='Paleomag Step') #Polarities, 0 false, 1 true
        #print 'mag_function', mag_function
        
        cdf = np.cumsum(mag_function)

        #''
#        print cdf[-1]
#        if cdf[-1] == 0:
#          plt.title('%s %s\nmu: %0.2f sd: %0.2f' % (currentevent, Polarity[currentevent], mu, sigma))
#          plt.ylim(-1, 2)
#          plt.show()
#          plt.close()
        #''

        if cdf[-1] == 0:
          print '\n', currentevent
          print 'range_strat', range_strat
          print 'mag_function', mag_function
          print 'pdf', pdf

        cdf /= np.max(cdf) # Normalize to one

        #plt.plot(range_strat, cdf, lw=2, label='CDF')
        #plt.title('From mag only, no radiometric age')
        #plt.show()
        #plt.close()

        SampledAges[currentevent] = np.min(np.where( cdf >= np.random.uniform() )) * Rdt + minage_strat
        #plt.axvline(SampledAges[currentevent], color='black', lw=2, label='Sampled Age')
        #plt.legend(loc='best')
        #plt.title('%s %s\nmu: %0.2f sd: %0.2f' % (currentevent, Polarity[currentevent], mu, sigma))
        #plt.ylim(-1, 2)
        #plt.show()
        #plt.close()


    if use_mag == True and sigma != -9999:
      #print currentevent, Polarity[currentevent]
      minage_strat = AcceptableAge_MIN
      maxage_strat = AcceptableAge_MAX
      #print 'Minage:', minage_strat, 'Maxage:', maxage_strat

      if Polarity[currentevent] == 'NA':
        minage = AcceptableAge_MIN
        maxage = AcceptableAge_MAX
        #print minage, maxage, mu, sigma
        # Convert minage and maxage to standard normal range because a, b are the standard deviations
        a = (minage_strat - mu) / sigma
        b = (maxage_strat - mu) / sigma
        #print a, b
        SampledAges[currentevent] = truncnorm.rvs(a, b, loc=mu, scale=sigma)
        if SampledAges[currentevent] <=0:
          while SampledAges[currentevent] <= 0:
            SampledAges[currentevent] = truncnorm.rvs(a, b, loc=mu, scale=sigma)

        if SampledAges[currentevent] == np.inf:
          #print 'Inf encountered at unit', currentevent
          SampledAges[currentevent] = np.random.uniform(maxage, minage)

      else:
        #x = np.linspace(-2.5, 2.5, 6)
        #import pylab as plt
        #np.random.seed(1234)
        range_strat = np.arange(minage_strat, maxage_strat, Rdt)
        #print 'Using this range', range_strat
        if Polarity[currentevent] == 'n': # Has radiometric age, polarity is 'normal'
          mag_function = mag_interpolate_n(range_strat)
          #plt.plot(range_strat, mag_function, label='Paleomagnetic Step Function') #Polarities, 0 false, 1 true
          #plt.title('Example of Paleomagnetic Step Function')
          #plt.ylabel('Probability')
          #plt.xlabel('Time (Ma)')
          #plt.ylim(-0.5, 1.5)
          #plt.xlim(range_strat[0], range_strat[-1])
          #plt.gca().invert_xaxis()
          #plt.show()
          #quit()

          #print 'mag_function', mag_function
          #print mu, sigma
          pdf = norm.pdf(range_strat, loc=mu, scale=sigma)
          #print pdf.max()*sigma
          #plt.plot(range_strat, pdf*sigma,label='PDF of Radiometric Age')
          #plt.title('Paleomag Step Function and PDF')
          #plt.gca().invert_xaxis()
          #plt.legend()
          #plt.show()
          #quit()
          #print 'pdf', pdf

          pdf_mag = pdf*mag_function
          #plt.plot(range_strat, pdf_mag*sigma)
          #plt.title('Gaussian PDF Times Paleomag Step Function')
          #plt.xlim(range_strat[0], range_strat[-1])
          #plt.gca().invert_xaxis()
          #plt.show()
          #quit()

          cdf = np.cumsum(pdf_mag)
          #print 'cdf', cdf
          #''
#          print cdf[-1]
#          if cdf[-1] == 0:
#            plt.title('%s %s\nmu: %0.2f sd: %0.2f' % (currentevent, Polarity[currentevent], mu, sigma))
#            plt.ylim(-1, 2)
#            plt.show()
#            plt.close()
          #''

          if cdf[-1] == 0:
            print '\n', currentevent
            print 'range_strat', range_strat
            print 'mag_function', mag_function
            print 'pdf', pdf

          cdf /= np.max(cdf) # Normalize to one

          #plt.plot(range_strat, cdf, label='CDF')
          #plt.xlim(range_strat[0], range_strat[-1])
          #plt.title('Normalized Cumulative Distribution Function')
          #plt.gca().invert_xaxis()
          #plt.show()
          #quit()
          #plt.close()
          
          #print np.min(np.where( cdf >= np.random.uniform() )) * Rdt + minage_strat
          rand_val = np.random.uniform()
          SampledAges[currentevent] = np.min(np.where( cdf >= rand_val )) * Rdt + minage_strat
          #plt.axhline(rand_val, color='red', lw=1, label='Random Value')
          #plt.axvline(SampledAges[currentevent], color='black', lw=1, label='Time Interval Sampled')
          #plt.legend()
          #plt.axvline(SampledAges[currentevent]-Rdt, color='black', lw=1)
          #plt.xlim(range_strat[0], range_strat[-1])
          #plt.gca().invert_xaxis()
          #plt.title('%s %s\nmu: %0.2f sd: %0.2f' % (currentevent, Polarity[currentevent], mu, sigma))
          #plt.ylim(-1, 2)
          #plt.show()
          #quit()
          #plt.close()

        if Polarity[currentevent] == 'r': # Has radiometric age, polarity is 'reversed'
          mag_function = mag_interpolate_r(range_strat)
          #plt.plot(range_strat, mag_function, lw=2, label='Paleomag Step') #Polarities, 0 false, 1 true
          #print 'mag_function', mag_function

          pdf = norm.pdf(range_strat, loc=mu, scale=sigma)
          #plt.plot(range_strat, pdf, lw=2, label='PDF', linestyle='dashed') #


          pdf_mag = pdf*mag_function
          #plt.plot(range_strat, pdf_mag, lw=2, label='PDF*Paleomag')

          cdf = np.cumsum(pdf_mag)

          #''
#          print cdf[-1]
#          if cdf[-1] == 0:
#            plt.title('%s %s\nmu: %0.2f sd: %0.2f' % (currentevent, Polarity[currentevent], mu, sigma))
#            plt.ylim(-1, 2)
#            plt.show()
#            plt.close()
          #''

          if cdf[-1] == 0:
            print '\n', currentevent
            print 'range_strat', range_strat
            print 'mag_function', mag_function
            print 'pdf', pdf

          cdf /= np.max(cdf) # Normalize to one

          #plt.plot(range_strat, cdf, lw=2, label='CDF')

          SampledAges[currentevent] = np.min(np.where( cdf >= np.random.uniform() )) * Rdt + minage_strat
          #plt.axvline(SampledAges[currentevent], color='black', lw=2, label='Sampled Age')
          #plt.legend(loc='best')
          #plt.title('%s %s\nmu: %0.2f sd: %0.2f' % (currentevent, Polarity[currentevent], mu, sigma))
          #plt.ylim(-1, 2)
          #plt.show()
          #plt.close()

    '''
    #if use_mag == False and sigma != 0:
    if sigma != 0:
      #sys.stdout.write('f')
      try:
        len(Ages[currentevent]) > 1 # If there is more than one age/uncertainty reported per event, randomly choose one.
        agechoice = np.random.randint(0, len(Ages[currentevent]))
        mu = Ages[currentevent][agechoice]
        sigma = Uncertainty[currentevent][agechoice]
      except:
        TypeError
        mu = Ages[currentevent]
        sigma = Uncertainty[currentevent]

      if mu == -9999:
        #sys.stdout.write('.')
        minage = AcceptableAge_MIN
        maxage = AcceptableAge_MAX
        SampledAges[currentevent] = np.random.uniform(maxage, minage)
        #sys.stdout.write('%0.3f' % SampledAges[currentevent])
        continue
        #print 'min, sampled, max', minage,SampledAges[currentevent], maxage

      else:
        #sys.stdout.write(',')
        minage = AcceptableAge_MIN
        maxage = AcceptableAge_MAX

        # Convert minage and maxage to standard normal range because a, b are the standard deviations
        a = (minage - mu) / sigma
        b = (maxage - mu) / sigma
        #print 'Here are a and b', a, b, 'for unit', currentevent

        # Use truncated normal distribution to sample age, make sure it is greater than zero
        SampledAges[currentevent] = truncnorm.rvs(a, b, loc=mu, scale=sigma)
        #sys.stdout.write("a=%0.3f b=%0.3f mu=%0.3f sig=%0.3f min=%0.3f max=%0.3f age=%0.3f" % (a, b, mu, sigma, minage, maxage, SampledAges[currentevent]))
        breakpt = 0 #break after 10
        if SampledAges[currentevent] <=0:
          while ((SampledAges[currentevent] <= 0) and (breakpt < 10)):
            SampledAges[currentevent] = truncnorm.rvs(a, b, loc=mu, scale=sigma)
            breakpt += 1
        if np.isinf(SampledAges[currentevent]):
          while ((np.isinf(SampledAges[currentevent])) and (breakpt < 10)):
            SampledAges[currentevent] = truncnorm.rvs(a, b, loc=mu, scale=sigma)
            breakpt += 1
        
        #If the sample is too far along the tail, just throw the age model out!
        # and choose from a random uniform.
        if breakpt == 10:
          SampledAges[currentevent] = np.random.uniform(maxage, minage)
        #sys.stdout.write('%0.3f' % SampledAges[currentevent])
        #sys.stdout.write('(%0.2f %0.2f)%0.1f-%d' % (minage, maxage, SampledAges[currentevent], breakpt))
        # Set the while counter to zero so that while loop can exit if tends towards infinite loop
        whilecountmin = 0
        whilecountmax = 0

        if SampledAges[currentevent] < minage:
          print '\n\nEntering the minage while loop: (a-b) = %f' % (a-b) 
          while SampledAges[currentevent] < minage and whilecountmin < 100:
            SampledAges[currentevent] = truncnorm.rvs(a, b, loc=mu, scale=sigma)
            print SampledAges[currentevent], 'should be greater than', minage, SampledAges[currentevent]>minage
            if SampledAges[currentevent] == np.inf:
              #print 'Inf encountered at unit', currentevent
              SampledAges[currentevent] = np.random.uniform(maxage, minage)
            if b-a < 0.1: # and whilecountmin == 100:
              SampledAges[currentevent] = np.random.uniform(maxage, minage)
              print 'Difference between b and a less than 0.1:', b-a, currentevent
              print 'minage:%f sampledage:%f, maxage:%f' % (minage, SampledAges[currentevent], maxage)
              break
            whilecountmin += 1
          print SampledAges[currentevent], 'should be greater than', minage, SampledAges[currentevent]>minage

        if SampledAges[currentevent] > maxage:
          print '\n\nEntering the maxage while loop: (a-b) = %f' % (a-b)
          while SampledAges[currentevent] > maxage and whilecountmax < 100:
            SampledAges[currentevent] = truncnorm.rvs(a, b, loc=mu, scale=sigma)
            #print SampledAges[currentevent], 'should be less than', maxage, SampledAges[currentevent]<maxage
            if SampledAges[currentevent] == np.inf:
              #print 'Inf encountered at unit', currentevent
              SampledAges[currentevent] = np.random.uniform(maxage, minage)
            if b-a < 0.1: #and whilecountmax == 100:
              SampledAges[currentevent] = np.random.uniform(maxage, minage)
              #print 'Difference between b and a less than 0.1:', b-a, currentevent
              #print 'minage:%f sampledage:%f, maxage:%f' % (minage, SampledAges[currentevent], maxage)
              break
            whilecountmax += 1
          #print SampledAges[currentevent], 'should be less than', maxage, SampledAges[currentevent]<maxage
        
        '''
        if whilecountmin != 0:
          print 'It took %i tries within the minage while loop for event: %s' % (whilecountmin, currentevent)

        if whilecountmax != 0:
          print 'It took %i tries within the maxage while loop for event: %s' % (whilecountmax, currentevent)
        '''

        # Can optimize by removing these and only generating them by reruning this run idx if necessary
        statement_min.append('%f > %f %s' % (SampledAges[currentevent], minage, SampledAges[currentevent]>minage))
        statement_min.append('a:%f b:%f (a-b):%f' % (a, b, a-b))
        statement_min.append('(%f - %f) / %f' % (minage, mu, sigma))
        statement_min.append('minage:%f maxage:%f max-min:%f' % (minage, maxage, maxage-minage))
        statement_max.append('%f < %f %s' % (SampledAges[currentevent], maxage, SampledAges[currentevent]<maxage))
        statement_max.append('a:%f b:%f (a-b):%f' % (a, b, a-b))
        statement_max.append('(%f - %f) / %f' % (minage, mu, sigma))
        statement_max.append('minage:%f maxage:%f max-min:%f' % (minage, maxage, maxage-minage))

        if whilecountmin == 100 or whilecountmax == 100:
          print 'Aborting after 100 tries'
          print pastevents
          for entry in statement_min:
            print entry
          for entry in statement_max:
            print entry
          print 'Creating graph'
          SampledAgeskeys = np.array(SampledAges.keys())
          SampledAgesvalues = np.array(SampledAges.values())
          error_at = currentevent
          error_min = SampledAgeskeys[np.where(SampledAgesvalues == minage)]
          error_max = SampledAgeskeys[np.where(SampledAgesvalues == maxage)]
          #MakeStratGraph(runID, SampledAges, error_at, error_min, error_max, pastevents)      
          #quit()

        if SampledAges[currentevent] < minage:
          print 'This cannot be a minage issue!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          quit()

        if SampledAges[currentevent] > maxage:
          print 'This cannot be a maxage issue!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          quit()

  return SampledAges

def sample_ages1(events, relationships, runID, Ages, Uncertainty, event_ageKey=None):
	'''
	This function generates a dictionary of sampled ages. 
	It iterates over a list of events. 
	It finds an acceptable age range based upon stratigraphic relationships and SampledAges. 
	This acceptable age range then truncates the probability distribution function. 
	A cumulative distribution function is created so that a random number from 0-1 
		can be selected and the associated age sampled, from the cumulative distribution function. 
	In cases where minage/maxage from stratigraphy truncate the probability distribution function 
		to one of the tails where the probability of an event is zero, the code then chooses a 
		random uniform age between the two bounding units. 
	In this way, stratigraphy can override the radiometric date.
	'''
	SampledAges = {} # Need a for loop to create a new dictionary for SampledAges
	#print 'relationships', relationships
	for i in Ages:
		SampledAges[i] = -9999 # Initially set SampledAges to -9999
	pastevents = []
	statement_min = []
	statement_max = []
	
	for currentevent in events:
		
		#####FIND EVENT AGE RANGE######
		AcceptableAge_MIN, statement_min = minage_finder_debug(currentevent,relationships,SampledAges,statement_min)
		AcceptableAge_MAX, statement_max = maxage_finder_debug(currentevent,relationships,SampledAges,statement_max)
		
		#Test for valid age range
		#if AcceptableAge_MIN >= AcceptableAge_MAX:
			#bad_range = str(AcceptableAge_MAX-AcceptableAge_MIN)
			#print "\n	ERROR: NO ACCEPTABLE AGE RANGE FOR EVENT ("+bad_range+" yrs)"+currentevent
			#print AcceptableAge_MIN, AcceptableAge_MAX
		
		# This is where the fun begins
		if event_ageKey == None: # If we didn't use key_stratigraphic_unit for event sorting
			#sys.stdout.write('1')
			try:
				len(Ages[currentevent]) > 1 # If there is more than one age/uncertainty reported per event, randomly choose one.
				agechoice = np.random.randint(0, len(Ages[currentevent]))
				mu = Ages[currentevent][agechoice]
				sigma = Uncertainty[currentevent][agechoice]
				#sys.stdout.write('a')
			except TypeError:
				mu = Ages[currentevent]
				sigma = Uncertainty[currentevent]
				#sys.stdout.write('b')
		
		if sigma == 0: # If the date is historic or exact...
			#sys.stdout.write('3')
			SampledAges[currentevent] = mu
			continue

		if sigma != 0:
			#sys.stdout.write('f')
			try:
				len(Ages[currentevent]) > 1 # If there is more than one age/uncertainty reported per event, randomly choose one.
				agechoice = np.random.randint(0, len(Ages[currentevent]))
				mu = Ages[currentevent][agechoice]
				sigma = Uncertainty[currentevent][agechoice]
			except:
				TypeError
				mu = Ages[currentevent]
				sigma = Uncertainty[currentevent]

			if mu == -9999:
				#sys.stdout.write('.')
				minage = AcceptableAge_MIN
				maxage = AcceptableAge_MAX
				SampledAges[currentevent] = np.random.uniform(maxage, minage)
				#sys.stdout.write('%0.3f' % SampledAges[currentevent])
				continue
				#print 'min, sampled, max', minage,SampledAges[currentevent], maxage

			else:
				#sys.stdout.write(',')
				minage = AcceptableAge_MIN
				maxage = AcceptableAge_MAX

				# Convert minage and maxage to standard normal range because a, b are the standard deviations
				a = (minage - mu) / sigma
				b = (maxage - mu) / sigma
				#print 'Here are a and b', a, b, 'for unit', currentevent
				
				#if b-a < 0.1
				#Just choose random normal if b-a < 0.1!!
				
				# Use truncated normal distribution to sample age, make sure it is greater than zero
				SampledAges[currentevent] = truncnorm.rvs(a, b, loc=mu, scale=sigma)
				#sys.stdout.write("a=%0.3f b=%0.3f mu=%0.3f sig=%0.3f min=%0.3f max=%0.3f age=%0.3f" % (a, b, mu, sigma, minage, maxage, SampledAges[currentevent]))
				breakpt = 0 #break after 10
				if SampledAges[currentevent] <=0:
					while ((SampledAges[currentevent] <= 0) and (breakpt < 10)):
						SampledAges[currentevent] = truncnorm.rvs(a, b, loc=mu, scale=sigma)
						breakpt += 1
				if np.isinf(SampledAges[currentevent]):
					while ((np.isinf(SampledAges[currentevent])) and (breakpt < 10)):
						SampledAges[currentevent] = truncnorm.rvs(a, b, loc=mu, scale=sigma)
						breakpt += 1
				
				#If the sample is too far along the tail, just throw the age model out!
				# and choose from a random uniform.
				if breakpt == 10:
					SampledAges[currentevent] = np.random.uniform(maxage, minage)

				#Append 
				if SampledAges[currentevent] < minage:
					print 'This might be a minage issue!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
					return -1

				if SampledAges[currentevent] > maxage:
					print 'This might be a maxage issue!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
					return -1

	return 0

###############################################################################
def getages(Ages, style, numruns, relationships, Uncertainty):
	'''
	getages sorts events based on user-input strategies, and runs numruns # of 
	simulations to sample event ages with sample_ages. Timing information is kept
	of sorting and sampling elapsed time
	'''
	
	#Result, timing, and list variables that are pre-allocated
	results = np.zeros((len(Ages), numruns))# space to save results
	SampledAgeslist = []      # lists of keys from SampledAges to check sorting of results
	eventtiming = np.zeros(numruns) #Time it takes for sorting events
	resulttiming = np.zeros(numruns) #Time it takes for sampling ages
	
	### Create Pre-sorting lists for some methods
	if style == 'topdown':
		tdlist = make_topdown_list()
	elif style == 'bottomup':
		bulist = make_bottomup_list()
	elif style == 'outside_in':
		top1 = make_topdown_list()[0]
		bottom1 = make_bottomup_list()[0]
	elif style == 'user_defined':
		print 'Events will be randomly shuffled within each list, from top to bottom'
		for idx,sublist in enumerate(Order_list):
			print 'Level', idx+1, sublist
	
	### VEAM Simulations begin here:
	for i in range(numruns):
		eventstart = time.time()
		event_ageKey = None
		
		### Determine Sorting Method
		if style == 'random':
			#Takes about 0.72 seconds per run for 10 and 100 runs
			events = random_events(Ages)
		elif style == 'topdown':
			#Takes about 0.97 seconds per run for 100 runs. 
			#Perhaps optimize the events list algorithm by saving the stratigraphic 
			#orders as a list of lists, then random shuffle each list within the list 
			#and assemble new eventlist? Not really a lot of work to spare, 
			#could take more time...
			events = topdown_events(tdlist)
		elif style == 'bottomup':
			#Takes about 0.61 seconds per run for 100 runs
			events = bottomup_events(bulist)
		elif style == 'outside_in':
			#This takes about #### seconds per run for 100 runs
			events = outside_in_events(top1, bottom1)
		elif style == 'most_contacts':
			#This takes about 0.85 seconds per run for 100 runs
			events = most_contacts_events(Ages, relationships)
		elif style == 'most_contacts_list':
			#This takes about 0.77 seconds per run for 100 runs. 
			#It has more code than most_contacts... 
			#Must have to do with sample_ages and the way the recursive function works???
			events = most_contacts_list_events()
		elif style == 'crater_age_uncertainty':
			#This takes about 0.96 seconds per run for 100 runs
			events = crater_age_uncertainty_events(Uncertainty)
		elif style == 'user_defined':
			events = user_defined_events()
		elif style == 'key_stratigraphic_unit':
			#This takes about 0.96 seconds per run for 100 runs
			events, event_ageKey = key_stratigraphic_unit_events()
		elif style == 'ignore_strat':
			#This takes about 1.64 seconds per run for 100 runs
			relationships = [[0, 0], [0, 0]]
			events = random_events(Ages)
		eventtiming[i] = time.time() - eventstart # Time elapsed to create sorted event list
		
		### Sample Ages
		resultstart = time.time()
		sys.stdout.write("\rOn run %d out of %d" % ((i+1),numruns)) # i+1, 'runs out of', numruns
		result = sample_ages(events, relationships, i, Ages, Uncertainty, event_ageKey)
		resulttiming[i] = time.time() - resultstart
		
		#Store results of events and their sampled ages
		for idx,k in enumerate(result):
			#print i, idx, k, result[k]
			#if min(result) != -9999:
			results[idx, i] = result[k]
		SampledAgeslist.append(result.keys())
	
	# Error checking
	for i in range(1, len(SampledAgeslist)):
		if SampledAgeslist[i-1] != SampledAgeslist[i]:
			print '\nYa done goofed up... Not all data are in the proper order\n'
			quit()
	
	#Done with VEAM Simulations
	end = time.time()
	print '\n\nThis %s simulation took' % (style), end - start, 'seconds\n'
	print 'The event list creation took', np.mean(eventtiming), 'seconds'
	print 'The age sampling took', np.mean(resulttiming), 'seconds'

	#Take the transpose because infile data are referenced as results[0] = all ages for Event[0]
	resultsT = np.transpose(results)
	#print resultsT
	
	### Check for bad runs
	#Save the row numbers for runs that don't have -9999 in them#
	#Rows with -9999 had an error#
	nonzeroruns = []
	badruns = []
	for idx, i in enumerate(resultsT):
		if -9999 not in i:
			nonzeroruns.append(idx)
		if -9999 in i:
			badruns.append(idx)
			print 'this bad run idx', idx
			for j, ev in enumerate(i):
				#print '%s, %0.2f;' % (SampledAgeslist[0][j], ev),
				if ev < 0:
					print '%s, %0.2f\t' % (SampledAgeslist[0][j], ev),
			print '\ntotal badruns', badruns

	return results, badruns, SampledAgeslist[0]

def veam_main(inputs):
	#Inputs needs to be an instance of a Variable Class
	
	#Define global minimum/maximum model age
	
	global GLOBAL_MINAGE
	global GLOBAL_MAXAGE
	global start
	
	start = time.time()
	
	GLOBAL_MINAGE = inputs.minAge
	GLOBAL_MAXAGE = inputs.maxAge
	#Define other variables
	ages_db_file = inputs.ageDB
	strat_db_file = inputs.stratDB
	use_mag = inputs.geoMag
	Rdt = inputs.res
	numruns = inputs.sims
	style = inputs.sorting
	
	print '\nuse_mag:', use_mag
	
	field = eventLibrary()
	
	#cwd = os.getcwd()
	for i in range(0,1):
	  relationships, Ages, Uncertainty, Polarity, Order_list = load_databases(strat_db_file,ages_db_file)
	  
	  '''
	  Changing over to class!!
	  '''
	  for item in Ages:
				field.addEvent(item)
				field.events[-1].addAgeModel(Ages[item], Uncertainty[item]) #Does this add multiple ages? It must!!
				
	  field.sortAgeUncertainty()
	  
	  '''Done with this'''
	  
	  if use_mag == True:
	    Rdt = Rdt # This is the timestep used for the cumulative density function
	    pol_n, pol_r, tops, bases = load_mag_timescale()
	    #print tops[pol_n]
	    #print bases[pol_r]
	    #print pol_n
	    #print pol_r
	    #print tops
	    #print bases
	    # This section of code creates a 1d interpolation of the geomagnetic timescale, 0 = Normal, 1 = Reversed
	    xs = []
	    yr = []
	    yn = []
	    for topn, basen, topr, baser in zip(tops[pol_n], bases[pol_n], tops[pol_r], bases[pol_r]):
	      xs.append(topn)
	      xs.append(basen)
	      xs.append(topr)
	      xs.append(baser)
	
	      yr.append(0)
	      yr.append(0)
	      yr.append(1)
	      yr.append(1)
	
	      yn.append(1)
	      yn.append(1)
	      yn.append(0)
	      yn.append(0)
	
	    #print xs
	    #print yr
	    mag_interpolate_r = interp1d(xs, yr)
	    mag_interpolate_n = interp1d(xs, yn)
	    #xspl = np.arange(xs[0], xs[-1], Rdt)
	    #import pylab as plt
	    #plt.plot(xspl, mag_interpolate_r(xspl), lw=2, color='red')
	    #plt.plot(xspl, mag_interpolate_n(xspl), lw=2, color='green')
	    #plt.axvline(0.5)
	    #print mag_interpolate(0.5)
	    #print xspl
	    #stepfn = np.zeros(len(xspl))
	    #plt.plot(xs, ys)
	    #plt.ylim(-1, 2)
	    #plt.show()
	  #print Polarity
	  #quit()
	  '''
	  print '\nHere are the stratigraphic relationships\n', relationships
	  print '\nHere are the ages\n', Ages
	  print '\nHere are the uncertainties\n',Uncertainty
	  '''
	  #quit()
	  print '######################\nWorking on style', style, numruns
	  results, badruns, SampledAgeslist = getages(Ages=Ages, style=style, numruns=numruns, relationships=relationships, Uncertainty=Uncertainty) # More than just the badruns comes back
	  
	  resultsT = np.transpose(results)
	
	  ####Save the row numbers for runs that don't have -9999 in them####
	  ####Rows with -9999 had an error####
	  nonzeroruns = []
	  badruns = []
	  for idx, resultT in enumerate(resultsT):
	    if -9999 not in resultT:
	      nonzeroruns.append(idx)
	    if -9999 in resultT:
	      badruns.append(idx)
	      #print 'resultT', resultT 
	      #print 'idx', idx
	      #print 'badruns', badruns
	
	
	  ####These are the end results####
	  endresultsT = resultsT[nonzeroruns] 
	  endresults = np.transpose(endresultsT) 
	  print 'There are', len(nonzeroruns), 'successful runs out of', np.shape(results)[1]
	 
	  # Save the results
	  np.save('results_%s' % (style), endresults)
	
	  if len(nonzeroruns) == 0:
	    print 'No successful runs'
	    return -1
	
	np.save('%s_eventlist.npy' % (style), SampledAgeslist)
	
	datafile = 'results_%s.npy' % (style)
	data = np.load(datafile)
	data = np.transpose(data)
	
	Eventlist = np.load('%s_eventlist.npy' % (style))
	
	print "Eventlist",Eventlist
	
	#if len(nonzeroruns) > 0:
	#print len(data[0])
	
	
	out = open('%s_ASCII.txt' % style, 'w')
	print >> out, 'Parameters Used To Generate The Data Herein'
	print >> out, 'Age Database: %s/%s' % (os.getcwd(), ages_db_file)
	print >> out, 'Strat Database: %s/%s' % (os.getcwd(), strat_db_file)
	print >> out, 'S: %0.5f' % GLOBAL_MAXAGE
	print >> out, 'T: %0.5f' % GLOBAL_MINAGE
	print >> out, 'Rdt: %0.5f' % Rdt
	print >> out, 'Use Mag: %s' % use_mag
	for idx, event in enumerate(Eventlist):
	  if idx != len(Eventlist)-1:
	    print >> out, '%s,' % event,
	  else:
	    print >> out, event
	
	for run in range(np.shape(data)[0]):
	  #print run
	  for idx, age in enumerate(data[run]):
	    if idx != len(Eventlist)-1:
	      print >> out, '%s,' % age,
	    else:
	      print >> out, age
	
	out.close()
	
	return 0

class event():
	def __init__(self):
		self.id = ''
		self.ageModel    = [] # list of modeled ages
		self.uncertModel = [] # list of associated uncertainties
		self.veamAges    = [] # VEAM-simulated ages
		self.polarity    = None #Normal='N', Reverse='R'?
		self.stratAbove  = [] # ids of events that are immediately stratigraphically above
		self.stratBelow  = [] # ids of events that are immediately stratigraphically lower
		self.totalRels   = 0  # total number of stratigraphic relationships
		self.modelChoice = 0  #
		self.sortOrder   = 0
		
	def addAgeModel(self, time, uncertainty):
		self.ageModel.append(float(time))
		self.uncertModel.append(float(uncertainty))
	
	def checkStratRels(self):
		for rel in self.stratBelow:
			if rel in self.stratAbove:
				return False #if a strat relationship is in both upper and lower lists, that's an error
		return True #these lists are ok.
		
		def calcTotalRels(self):
			self.totalRels = len(self.stratAbove) + len(self.stratBelow)
		
		def chooseAgeModel(self):
			if len(self.ageModel) > 1:
				self.modelChoice = np.random.randint(len(self.ageModel))
	
class eventLibrary():
	def __init__(self):
		self.events = []

	def addEvent(self,name):
		### Create a new event
		e = event()
		e.id = str(name)
		self.events.append(e)
		
	def checkStrat(self):
		### Check Stratigraphy Web. 
		#   If valid, return True, else return False
		for e in self.events:
			check = e.checkStratRels
			if check == False:
				return False
		
		#make a deeper check
		
		return True
		
	
	def sortRandom(self):
		### Sorts all events randomly
		np.random.shuffle(self.events)
		
		#Choose a random age model if event has multiple
		for e in self.events:
			e.chooseAgeModel()
	
	def sortMostContacts(self):
		### Sorts all events so events with more contacts are dated first.
		# Events with the same number of contacts get a random shuffle
		sortCount = 0 #Will be assigned in order 0 -> N for each event
		
		for e in self.events:
			#Choose a random age model if event has multiple
			e.chooseAgeModel()
			#Calculate number of total strat relationships for each event
			e.calcTotalRels()
		
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
			
		#Choose a random age model if event has multiple
		for e in self.events:
			e.chooseAgeModel()
	
	def sortAgeUncertainty(self):
		### Sorts all events so events with least age model uncertainty are dated first.
		# Events with the same age uncertainty will get a random shuffle.
		# Events with more than one age model will have one chosen for them.
		sortCount = 0 #Will be assigned in order 0 -> N for each event
		
		#Choose a random age model if event has multiple
		for e in self.events:
			e.chooseAgeModel()
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
		self.ageDB   = ''
		self.stratDB = ''
		self.res     = 0.0
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
		

class Window(QtWidgets.QWidget):
	
	def __init__(self):
		super(Window, self).__init__()
		self.init_ui()
		
	def init_ui(self):
		self.chkct = 0
		
		### The Banner Logo for VEAM
		banner = QtWidgets.QHBoxLayout()
		self.bannerLabel = QtWidgets.QLabel()
		self.bannerLabel.setPixmap(QtGui.QPixmap('greeley_cone_banner_300dpi.png').scaled(600,
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
		self.leYoung  = QtWidgets.QLineEdit()
		self.leOld    = QtWidgets.QLineEdit()
		self.leRes    = QtWidgets.QLineEdit()
		self.cbxYoung = QtWidgets.QComboBox(self)
		self.cbxOld   = QtWidgets.QComboBox(self)
		self.cbxRes   = QtWidgets.QComboBox(self)
		#self.optYoung = QtWidgets.QComboBox
		
		self.cbxYoung.addItems(('Ma', 'ka', 'a'))
		self.cbxOld.addItems(('Ma', 'ka', 'a'))
		self.cbxRes.addItems(('Ma', 'ka', 'a'))
		
		hYoung = QtWidgets.QHBoxLayout() #Young Hor Entry Box
		hYoung.addWidget(self.leYoung)
		hYoung.addWidget(self.cbxYoung)
		hOld = QtWidgets.QHBoxLayout() #Old Hor Entry Box
		hOld.addWidget(self.leOld)
		hOld.addWidget(self.cbxOld)
		gAges = QtWidgets.QGridLayout() #Young-Old Grid
		gAges.setVerticalSpacing(2)
		gAges.setHorizontalSpacing(100)
		gAges.addWidget(self.lblYoung, 1, 1)
		gAges.addWidget(self.lblOld, 1, 2)
		gAges.addLayout(hYoung, 2, 1)
		gAges.addLayout(hOld, 2, 2)
		hRes = QtWidgets.QHBoxLayout() #Age Resolution Hor Box
		hRes.addWidget(self.lblRes)
		hRes.addWidget(self.leRes)
		hRes.addWidget(self.cbxRes)
		hRes.addStretch()
		
		vAges = QtWidgets.QVBoxLayout() #Total Age Range Box Layout
		vAges.addLayout(gAges)
		vAges.addLayout(hRes)
		
		gpAges = QtWidgets.QGroupBox('Potential Volcanic Age Range') #Age Span groupbox
		gpAges.setLayout(vAges)
		
		### File Paths
		self.lblAgeDB   = QtWidgets.QLabel('Age Database Path ')
		self.lblStratDB = QtWidgets.QLabel('Stratigraphy Database Path ')
		self.leAgeDB    = QtWidgets.QLineEdit()
		self.leStratDB  = QtWidgets.QLineEdit()
		self.bAgeDB     = QtWidgets.QPushButton('...')
		self.bStratDB   = QtWidgets.QPushButton('...')
		self.cbxAgeUnits = QtWidgets.QComboBox()
		
		self.cbxAgeUnits.addItems(('Ma', 'ka', 'a'))
		
		gPath = QtWidgets.QGridLayout() #File Paths group Layout
		gPath.setSpacing(2)
		gPath.addWidget(self.lblAgeDB, 1, 0)
		gPath.addWidget(self.leAgeDB, 1, 1)
		gPath.addWidget(self.bAgeDB, 1, 2)
		gPath.addWidget(self.cbxAgeUnits, 1, 3)
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
		self.setGeometry(100, 100, 600, 400)
		
		### Actions
		
		self.bSaveCfg.clicked.connect(self.saveCfg)
		self.bLoadCfg.clicked.connect(self.loadCfg)
		self.bCheck.clicked.connect(lambda: self.checkCfg())
		self.bSubmit.clicked.connect(self.runVeam)
		
		self.bAgeDB.clicked.connect(lambda: self.open_path(self.leAgeDB))
		self.bStratDB.clicked.connect(lambda: self.open_path(self.leStratDB))
		
		self.show()
	
	def runVeam(self):
		### Run a check function that formats the form data first!
		
		
		check = 0
		check = self.checkCfg()
		
		if check == 0:
			veamVars = variables()
			veamVars.setMaxAge(self.vMaxAge)
			veamVars.setMinAge(self.vMinAge)
			veamVars.setResolution(self.vResol)
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
			self.bSubmit.setStyleSheet('QPushButton {background-color: red; font-weight: bold;}')
			self.bSubmit.repaint()
			
			ret = veam_main(veamVars)
			if ret==0:
				self.lblSubmit.setText('VEAM Completed!!')
				print 'VEAM Completed!!'
			else:
				self.lblSubmit.setText('VEAM had an error!')
			self.bSubmit.setStyleSheet('QPushButton {background-color: skyblue; font-weight: bold;}')
		
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
		
		errStyle = 'QLineEdit {background-color: pink}'
		valStyle = 'QLineEdit {background-color: white}'
		errMsg   = []
		
		#Find ComboBox Units (Ma, ka, or a)
		unitStrYoung = str(self.cbxYoung.currentText())
		unitStrOld   = str(self.cbxOld.currentText())
		unitStrRes   = str(self.cbxRes.currentText())
		self.unitStrAgeDB = str(self.cbxAgeUnits.currentText())
		if self.unitStrAgeDB=='Ma': multDB = 1e-6;
		elif self.unitStrAgeDB=='ka': multDB = 1e-3;
		else: multDB = 1; #unit will be 'a'
		if unitStrYoung=='Ma': multYoung = 1e6 * multDB;
		elif unitStrYoung=='ka': multYoung = 1e3 * multDB;
		else: multYoung = multDB;
		if unitStrOld=='Ma': multOld = 1e6 * multDB;
		elif unitStrOld=='ka': multOld = 1e3 * multDB;
		else: multOld = multDB;
		if unitStrRes=='Ma': multRes = 1e6 * multDB;
		elif unitStrRes=='ka': multRes = 1e3 * multDB;
		else: multRes = multDB;
		
		try:                                  #Young lineEdit
			minA = float(self.leYoung.text()) * multYoung #Call variable to float
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
			maxA = float(self.leOld.text()) * multOld
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
			res = float(self.leRes.text()) * multRes
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
			self.cbxYoung.setCurrentText(self.unitStrAgeDB)
			self.cbxOld.setCurrentText(self.unitStrAgeDB)
			self.cbxRes.setCurrentText(self.unitStrAgeDB)
		except IndexError: #If no unit is found, just assume 'Ma'
			self.cbxAgeUnits.setCurrentText('Ma')
			self.cbxYoung.setCurrentText('Ma')
			self.cbxOld.setCurrentText('Ma')
			self.cbxRes.setCurrentText('Ma')
		
		self.checkCfg()
		
		self.lblWarn.setText('Configuration Loaded!')


def main():
	### Define QT App Instance
	app = QtCore.QCoreApplication.instance()
	if app is None:
		app = QtWidgets.QApplication(sys.argv)

	### Create a Window
	VEAM_window = Window()
	
	#Exit upon closed window
	sys.exit(app.exec_())

if __name__ == '__main__':
	main()
