#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 16:42:04 2017
@author: James Wilson, Jacob Richardson

VEAM is a VOLCANIC EVENT AGE MODELER
"""

import sys
from PyQt5 import QtWidgets, QtCore, QtGui
from numpy import genfromtxt

### VEAM Variables and Structures ###

class Event():
	def __init__(self):
		self.id = ''
		self.stratAbove  = [] # ids of events that are immediately stratigraphically above
		self.stratBelow  = [] # ids of events that are immediately stratigraphically lower
		self.allAbove    = [] # ids of events indirectly above this event
		self.allBelow    = [] # ids of events indirectly below this event
		self.allAboveInd = [] # event library indices of allAbove
		self.allBelowInd = [] # event library indices of allBelow
	
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
			eID = e.id
			#print("Interrogating event:", eID)
			errorFlag = self.findProblemStrat(eID, relList, [eID]) 
			if len(errorFlag):
				return errorFlag
		return []
		
		#make a deeper check
		
		return True
	
	def findProblemStrat(self,curEvent,relList,eventsAbove):
		#finds problematic Contridictions in the Stratigraphy Web.
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

### GUI WINDOW ###
class mainWindow(QtWidgets.QWidget):
	
	def __init__(self, parent=None):
		super(mainWindow, self).__init__()
		self.init_ui()
		
	def init_ui(self):
		self.eventLib = eventLibrary()
		self.stratErrorList = []
		self.relList = []
		self.justSaved = True
		
		### The Banner Logo for VEAM
		banner = QtWidgets.QHBoxLayout()
		self.bannerLabel = QtWidgets.QLabel()
		self.bannerLabel.setPixmap(QtGui.QPixmap('richardson_SP_banner_600px.png').scaled(600,
												600, QtCore.Qt.KeepAspectRatio))
		banner.addStretch()
		banner.addWidget(self.bannerLabel)
		banner.addStretch()
		
		self.lblGreeley = QtWidgets.QLabel('<i>Tool Design: Jacob Richardson, 2020</i>')
		hCredit = QtWidgets.QHBoxLayout() #Young-Old Hor Box
		hCredit.addStretch()
		hCredit.addWidget(self.lblGreeley)
		
		### Load and Save Buttons
		self.bLoadAges  = QtWidgets.QPushButton('Load Events Database')
		self.bLoadStrat = QtWidgets.QPushButton('Load Stratigraphy Database')
		self.bSaveStrat = QtWidgets.QPushButton('Save Stratigraphy Database')
		self.bClearRels = QtWidgets.QPushButton('Clear All Relationships')
		self.bClearAll = QtWidgets.QPushButton('Clear Everything')
		self.bClearAll.setStyleSheet('QPushButton {color: red;}')
		self.bClearRels.setDisabled(True)
		self.bClearAll.setDisabled(True)
		
		hLoadSave = QtWidgets.QHBoxLayout() #sims, geomag, strat use box
		hLoadSave.addWidget(self.bLoadAges)
		hLoadSave.addStretch()
		hLoadSave.addWidget(self.bLoadStrat)
		hLoadSave.addStretch()
		hLoadSave.addWidget(self.bSaveStrat)
		
		hClear = QtWidgets.QHBoxLayout() #sims, geomag, strat use box
		hClear.addWidget(self.bClearRels)
		hClear.addWidget(self.bClearAll)
		hClear.addStretch()
		
		### Assigned Column
		
		self.lblHigher = QtWidgets.QLabel('Higher Events')
		self.lblHigherCt = QtWidgets.QLabel('')
		self.listHigher = QtWidgets.QListWidget()
		self.lblLower = QtWidgets.QLabel('Lower Events')
		self.lblLowerCt = QtWidgets.QLabel('')
		self.listLower = QtWidgets.QListWidget()
		self.cbxCurEvent = QtWidgets.QComboBox()
		self.cbxCurEvent.setInsertPolicy(QtWidgets.QComboBox.InsertAlphabetically) 
		self.oldindex = 0 #Keeps track of previous event of interest
		
		hHigherLbls = QtWidgets.QHBoxLayout()
		hHigherLbls.addWidget(self.lblHigher)
		hHigherLbls.addStretch()
		hHigherLbls.addWidget(self.lblHigherCt)
		
		hLowerLbls = QtWidgets.QHBoxLayout()
		hLowerLbls.addWidget(self.lblLower)
		hLowerLbls.addStretch()
		hLowerLbls.addWidget(self.lblLowerCt)
		
		#entries = ['1','1','1','1','1','1','1']
		#self.listHigher.addItems(entries)
		self.listHigher.setVerticalScrollBar(QtWidgets.QScrollBar())
		self.listHigher.setSelectionMode(QtWidgets.QListWidget.ExtendedSelection)
		self.listHigher.setStyleSheet('QListWidget {background-color: #eee;}')
		
		self.listLower.setVerticalScrollBar(QtWidgets.QScrollBar())
		self.listLower.setSelectionMode(QtWidgets.QListWidget.ExtendedSelection)
		self.listLower.setStyleSheet('QListWidget {background-color: #eee;}')
		
		self.cbxCurEvent.addItem('Choose an event...')
		
		
		vAssigned = QtWidgets.QVBoxLayout()
		vAssigned.addLayout(hHigherLbls)
		vAssigned.addWidget(self.listHigher)
		vAssigned.addWidget(self.cbxCurEvent)
		vAssigned.addLayout(hLowerLbls)
		vAssigned.addWidget(self.listLower)
		
		gpAssigned = QtWidgets.QGroupBox('Stratigraphically Related Events') #Age Span groupbox
		gpAssigned.setLayout(vAssigned)
		
		
		### Arrow Button Column
		self.bAddHigher  = QtWidgets.QPushButton('<-')
		self.bAddHigher.setToolTip('Add event to higher list')
		self.bRemoveHigher  = QtWidgets.QPushButton('->')
		self.bRemoveHigher.setToolTip('Remove event from higher list')
		self.bAddLower  = QtWidgets.QPushButton('<-')
		self.bAddLower.setToolTip('Add event to lower list')
		self.bRemoveLower  = QtWidgets.QPushButton('->')
		self.bRemoveLower.setToolTip('Remove event from lower list')
		
		vArrowButtons = QtWidgets.QVBoxLayout()
		vArrowButtons.addStretch()
		vArrowButtons.addWidget(self.bAddHigher)
		vArrowButtons.addWidget(self.bRemoveHigher)
		vArrowButtons.addStretch()
		vArrowButtons.addWidget(self.bAddLower)
		vArrowButtons.addWidget(self.bRemoveLower)
		vArrowButtons.addStretch()
		
		### Unassigned Column 
		self.lblUnassignedCt = QtWidgets.QLabel('')
		self.listUnassigned = QtWidgets.QListWidget()
		self.listUnassigned.setVerticalScrollBar(QtWidgets.QScrollBar())
		self.listUnassigned.setSelectionMode(QtWidgets.QListWidget.ExtendedSelection)
		self.listUnassigned.setStyleSheet('QListWidget {background-color: #eee;}')
		
		hUnassignedLbl = QtWidgets.QHBoxLayout()
		hUnassignedLbl.addStretch()
		hUnassignedLbl.addWidget(self.lblUnassignedCt)
		
		vUnassigned = QtWidgets.QVBoxLayout()
		vUnassigned.addLayout(hUnassignedLbl)
		vUnassigned.addWidget(self.listUnassigned)
		gpUnassigned = QtWidgets.QGroupBox('Unrelated Events') #Age Span groupbox
		gpUnassigned.setLayout(vUnassigned)
		
		### Combine Columns
		hWorkspace = QtWidgets.QHBoxLayout() #Includes both main columns
		hWorkspace.addWidget(gpAssigned)
		hWorkspace.addLayout(vArrowButtons)
		hWorkspace.addWidget(gpUnassigned)
		
		### Footer - Warning and Relationship Box
		self.lblWarn = QtWidgets.QLabel('')
		self.lblWarn.setStyleSheet('QLabel {font-weight: bold;}')
		self.cbxRelList = QtWidgets.QComboBox()
		self.cbxRelList.addItem('No Stratigrapic Relationships')
		
		hWarnings = QtWidgets.QHBoxLayout() #Warning Flags Layout
		hWarnings.addWidget(self.cbxRelList)
		hWarnings.addStretch()
		hWarnings.addWidget(self.lblWarn)
		
		### Assemble in vertical box
		v_box = QtWidgets.QVBoxLayout() #Vertical box for all main window
		v_box.addLayout(banner)
		v_box.addLayout(hLoadSave)
		v_box.addLayout(hClear)
		v_box.addLayout(hWorkspace)
		v_box.addLayout(hWarnings)
		#v_box.addStretch()
		v_box.addLayout(hCredit)
		
		
		self.setLayout(v_box)
		self.setWindowTitle('VEAM | Volcanic Event Age Model | Strat Tool')
		self.setGeometry(200, 100, 600, 400)
		
		### Main Button Actions
		self.bLoadAges.clicked.connect(self.loadEvents)
		self.bLoadStrat.clicked.connect(self.loadStratFile)
		self.bSaveStrat.clicked.connect(self.saveStratFile)
		self.bClearRels.clicked.connect(self.clearRels)
		self.bClearAll.clicked.connect(self.clearAll)
		self.cbxCurEvent.currentIndexChanged.connect(self.selectCurEvent)
		
		#Arrow Button Actions
		self.bAddHigher.clicked.connect(lambda: self.moveEvent('addHigher'))
		self.bRemoveHigher.clicked.connect(lambda: self.moveEvent('removeHigher'))
		self.bAddLower.clicked.connect(lambda: self.moveEvent('addLower'))
		self.bRemoveLower.clicked.connect(lambda: self.moveEvent('removeLower'))
		
		self.show()

	
	def loadEvents(self,filename='unassigned',cols=1):
		### Clear previous work
		self.cbxCurEvent.clear()
		self.cbxCurEvent.addItem('Choose an event...')
		self.listHigher.setStyleSheet('QListWidget {background-color: #eee;}')
		self.listLower.setStyleSheet('QListWidget {background-color: #eee;}')
		self.listUnassigned.setStyleSheet('QListWidget {background-color: white;}')
		self.listUnassigned.clear()
		self.listHigher.clear()
		self.listLower.clear()
		
		### Open File if filename not alread assigned
		if filename==False:
			filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Database File Path')
		
		#load ages database using genfromtxt
		#VEAM dates volcanic events using the age database, so all events will
		#be in this file already
		try:
			events = genfromtxt(filename[0], skip_header=1,delimiter=',', dtype='unicode')
		except ValueError:
			sys.stderr.write('\n ERROR: Check for extra spaces, commas in the names of events in %s and try again\n' % filename)
			return -1, None
		except OSError:
			return 0
		
		#Reset Event Library
		self.eventLib = eventLibrary()
		eventIDs = []
		
		# Add events to event library
		for i,entry in enumerate(events):
			for col in range(cols):
				exist = 0
				idname = entry[col].strip()
				
				if len(self.eventLib.events) > 0:
					for e in self.eventLib.events:
						if e.id==idname: #if event already exists, just add the new age model
							exist = 1
							break
				if exist == 0: #if this is a new event
					self.eventLib.addEvent(idname)
					eventIDs.append(idname)
				
		### Add units to unnassigned column
		self.listUnassigned.addItems(sorted(eventIDs,key=str.casefold))
		self.cbxCurEvent.addItems(sorted(eventIDs,key=str.casefold))
			
		self.updateListCts()
		self.makeRelList()
		self.lblWarn.setText(str(len(eventIDs))+' Events Loaded!')
		sys.stdout.write('Loaded Event List from file:\n  '+filename[0])
		sys.stdout.write('\nEvent List has '+str(len(eventIDs))+' events.\n\n')
		
		#Let a user clear everything and reset the Just Saved Flag since there should be no relationships
		self.bClearAll.setDisabled(False)
		self.justSaved = True		

	
	def loadStratFile(self):
		### Open File
		filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Database File Path')
		
		#If there aren't any events, then load these relationships as the events
		if len(self.eventLib.events) == 0:
			self.loadEvents(filename=filename,cols=2)
		
		#Load relationships from file
		try:
			relationships = genfromtxt(filename[0], skip_header=1,delimiter=',', dtype='unicode')
		except ValueError:
			sys.stderr.write('\n ERROR: Check for extra spaces, commas in the names of events in %s and try again\n' % filename)
			return -1, None
		except OSError:
			return 0
		
		#add all relationships.
		missingRels = 0
		for rel in relationships:
			lowerID = rel[0].strip()
			higherID = rel[1].strip()
			#Check that all relationships are events in the event tree
			eventFound = 0
			for e in self.eventLib.events:
				if e.id == lowerID:
					eventFound += 1
					next
				if e.id == higherID:
					eventFound += 2
					next
				if eventFound ==3:
					#Found both events! Add the relationship
					self.addNewRelationship(lowerID,higherID)
					break
			if eventFound != 3:
				#If there's a missing relationship, add it to the tally
				missingRels += 1
		
		#If there were missing relationships, show a warning dialog box
		if missingRels:
			missMessage = "Warning: Events from %d Relationship(s) were not found in the current event library." % missingRels
			missMessage = missMessage+" These relationships were ignored."
			self.warnBox = QtWidgets.QMessageBox()
			self.warnBox.setIcon(QtWidgets.QMessageBox.Warning)
			self.warnBox.setText(missMessage)
			self.warnBox.setWindowTitle("Load Warning")
			self.warnBox.setStandardButtons(QtWidgets.QMessageBox.Ok)
			
			self.warnBox.exec_()
			
			#Make Just Saved Flag False, since relationships are not the same as the loaded file
			self.justSaved = False
		else:
			#Make Just Saved Flag True, since relationships are all the same as the loaded file 
			self.justSaved = True
		sys.stderr.write('\nLoaded Relationships from file:\n  %s\n' % filename[0])
		
		self.makeRelList()
		self.checkRelList()
				
		
	def clearAll(self):
		if self.justSaved == False:
			#Check with warning box to see if user really wants this
			self.warnBoxCA = QtWidgets.QMessageBox()
			self.warnBoxCA.setIcon(QtWidgets.QMessageBox.Warning)
			self.warnBoxCA.setText("There are unsaved relationships. Do you want to continue?")
			self.warnBoxCA.setWindowTitle("Clear Warning")
			self.warnBoxCA.setStandardButtons(QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
			
			retval = self.warnBoxCA.exec_()
			if retval == QtWidgets.QMessageBox.No:
				return 0
		
		#CLEAR EVERYTHING
		self.eventLib = eventLibrary()
		self.stratErrorList = []
		self.relList = []
		self.selectCurEvent()
		self.justSaved = True
		self.cbxCurEvent.clear()
		self.makeRelList()
		self.lblWarn.setText('All Cleared')
		self.bClearAll.setDisabled(True)
			
	
	def clearRels(self):
		if self.justSaved == False:
			#Check with warning box to see if user really wants this
			self.warnBoxCR = QtWidgets.QMessageBox()
			self.warnBoxCR.setIcon(QtWidgets.QMessageBox.Warning)
			self.warnBoxCR.setText("There are unsaved relationships. Do you want to continue?")
			self.warnBoxCR.setWindowTitle("Clear Warning")
			self.warnBoxCR.setStandardButtons(QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
			
			retval = self.warnBoxCR.exec_()
			if retval == QtWidgets.QMessageBox.No:
				return 0
		
		#ok, you gave them a chance, now clear all those relationships!
		for rel in self.relList:
			self.removeRelationship(rel[0], rel[1])
		
		#reset boxes
		self.makeRelList()
		self.selectCurEvent()
		self.justSaved = True
		self.lblWarn.setText('Relationships Cleared')
		
		
	def saveStratFile(self):
		### Run a check function for the stratigraphy
		self.checkRelList()
		
		### Select a file and save
		filename = QtWidgets.QFileDialog.getSaveFileName(self, 'Save Configuration File')
		try:
			with open(filename[0], 'w') as cfg:
				cfg.write('Lower, Higher\n')
				for rel in self.relList:
					cfg.write(rel[0]+', '+rel[1]+'\n')
		except EnvironmentError: # parent of IOError, OSError *and* WindowsError where available
			self.lblWarn.setText('Relationships not saved')
			return -1
		
		#If there's still an error, print a warning at least
		if len(self.stratErrorList) > 0: #Print Result in Orange Warning Label
			self.lblWarn.setText('Warning: Errors found, Saved anyway')
			sys.stdout.write('\nRelationships (WITH ERRORS) saved to file:\n  '+filename[0]+'\n')
		else:
			self.lblWarn.setText('Relationships Saved!')
			sys.stdout.write('\nRelationships saved to file:\n  '+filename[0]+'\n')
			
		
		self.justSaved = True
		return 0


	def moveEvent(self,direction):
		### Remove Lower to Unassigned
		if direction=='removeLower':
			# Remove relationship(s) from event library
			for relEvent in self.listLower.selectedItems():
				lowerID = relEvent.text()
				higherID = self.cbxCurEvent.currentText()
				relError = self.removeRelationship(lowerID,higherID)
				if relError != 0: print('New Relationship Error')
			
			#move events to other side of GUI
			# sort selected rows in descending order in order to compensate shifting due to takeItem
			rows = sorted([index.row() for index in self.listLower.selectedIndexes()],
				 reverse=True)
			for row in rows:
				self.listUnassigned.addItem(self.listLower.takeItem(row))
			self.listUnassigned.sortItems()
				
		### Add Lower from Unassigned
		elif direction=='addLower':
			# Add new relationship(s) to event library
			for relEvent in self.listUnassigned.selectedItems():
				lowerID = relEvent.text()
				higherID = self.cbxCurEvent.currentText()
				relError = self.addNewRelationship(lowerID,higherID)
				if relError != 0: print('New Relationship Error')
			
			#move events to other side of GUI
			# sort selected rows in descending order 
			rows = sorted([index.row() for index in self.listUnassigned.selectedIndexes()],
				 reverse=True)
			for row in rows:
				#print(self.listUnassigned.itemText(row))
				self.listLower.addItem(self.listUnassigned.takeItem(row))
			self.listLower.sortItems()
		
		### Remove Higher to Unassigned	
		elif direction=='removeHigher':
			# Remove relationship(s) from event library
			for relEvent in self.listHigher.selectedItems():
				higherID = relEvent.text()
				lowerID = self.cbxCurEvent.currentText()
				relError = self.removeRelationship(lowerID,higherID)
				if relError != 0: print('New Relationship Error')
			
			#move events to other side of GUI
			# sort selected rows in descending order
			rows = sorted([index.row() for index in self.listHigher.selectedIndexes()],
				 reverse=True)
			for row in rows:
				self.listUnassigned.addItem(self.listHigher.takeItem(row))
			self.listUnassigned.sortItems()
				
		### Add Higher from Unassigned
		elif direction=='addHigher':
			# Add new relationship(s) to event library
			for relEvent in self.listUnassigned.selectedItems():
				higherID = relEvent.text()
				lowerID = self.cbxCurEvent.currentText()
				relError = self.addNewRelationship(lowerID,higherID)
				if relError != 0: print('New Relationship Error')
				
			# Move events to other side of GUI
			# sort rows in descending order 
			rows = sorted([index.row() for index in self.listUnassigned.selectedIndexes()],
				 reverse=True)
			for row in rows:
				self.listHigher.addItem(self.listUnassigned.takeItem(row))
			self.listHigher.sortItems()
			
		
		#Update Event Counts above each list box after every move
		self.updateListCts()
		relListError = self.makeRelList()
		if relListError != 0: print('Relationship List Error')
		
		#Look for infinite Loops
		self.checkRelList()
		
		#Change the Just Saved Flag
		self.justSaved = False
		
			
	def addNewRelationship(self,lowerID,higherID):
		#Add Stratigraphic Relationships to Events
		#event.stratAbove and stratBelow
		both = 0
		for event in self.eventLib.events:
			#if the event is the lower event, append the higher event to stratAbove
			if event.id == lowerID:
				event.stratAbove.append(higherID)
				both += 2
			#if the event is the higher event, append the higher event to stratBelow
			elif event.id == higherID:
				event.stratBelow.append(lowerID)
				both += 1
			#If both events have been found, go to next relationship
			if both == 3:
				return 0
			
		return -1

	
	def removeRelationship(self,lowerID,higherID):
		#Add Stratigraphic Relationships to Events
		#event.stratAbove and stratBelow
		both = 0
		for event in self.eventLib.events:
			#if the event is the lower event, append the higher event to stratAbove
			if event.id == lowerID:
				event.stratAbove.remove(higherID)
				both += 2
			#if the event is the higher event, append the higher event to stratBelow
			elif event.id == higherID:
				event.stratBelow.remove(lowerID)
				both += 1
			#If both events have been found, go to next relationship
			if both == 3:
				return 0
			
		return -1

	
	def makeRelList(self):
		#Make a 2xN array of all "Strat Below" relationships
		#Also check that "Strat Above" relationships are as numerous as stratbelow's
		stratsBelow = 0
		stratsAbove = 0 
		self.relList = []
		for event in self.eventLib.events:
			stratsAbove += len(event.stratAbove)
			stratsBelow += len(event.stratBelow)
			if len(event.stratBelow):
				for lowerID in event.stratBelow:
					self.relList.append([lowerID,event.id])
		if (stratsAbove-stratsBelow) != 0:
			return -1
		
		self.relList = sorted(self.relList)
		
		self.cbxRelList.clear()
		if len(self.relList) == 0:
			self.cbxRelList.addItem('No Stratigraphic Relationships')
			self.bClearRels.setDisabled(True)
		elif len(self.relList) == 1:
			self.cbxRelList.addItem('1 Stratigraphic Relationship')
			self.cbxRelList.addItem('Lower, Upper')
			self.cbxRelList.insertSeparator(2)
			self.cbxRelList.addItem(str(self.relList[0][0])+', '+str(self.relList[0][1]))
			
			#Let people click the clear relationship button
			self.bClearRels.setDisabled(False)
		else:
			self.cbxRelList.addItem(str(len(self.relList))+' Total Relationships')
			self.cbxRelList.addItem('Lower, Upper')
			self.cbxRelList.insertSeparator(2)
			for r in self.relList:
				self.cbxRelList.addItem(str(r[0])+', '+str(r[1]))
			
			#Let people click the clear relationship button
			self.bClearRels.setDisabled(False)
		
		self.cbxRelList.repaint()
		return 0
	
	
	def checkRelList(self):
		self.listUnassigned.setStyleSheet('QListWidget {background-color: white;}')
		self.listLower.setStyleSheet('QListWidget {background-color: white;}')
		self.listHigher.setStyleSheet('QListWidget {background-color: white;}')
		self.listUnassigned.repaint()
		self.listLower.repaint()
		self.listHigher.repaint()
		
		relListError = self.eventLib.checkAllStrat(self.relList)
		if len(relListError): 
			#Print out if the error is new
			if self.stratErrorList != relListError:
				sys.stderr.write("\n\nERROR: There is a stratigraphic loop with\n these vents:\n")
				for event in relListError:
					sys.stderr.write("  - "+event+"\n")
				self.lblWarn.setText('Relationship Loop Error!')
				self.stratErrorList = relListError
			
			#Paint everything white before things get painted magenta
			self.getRelsForCurEvent()
			
			#Find all bad events and paint magenta
			for badEvent in relListError:
				#Find bad entries in all lists
				bEIdx = self.listUnassigned.findItems(badEvent,QtCore.Qt.MatchExactly)
				if len(bEIdx) > 0:
					for B in bEIdx:
						row = self.listUnassigned.row(B)
						self.listUnassigned.item(row).setBackground(QtCore.Qt.magenta)
				bEIdx = self.listHigher.findItems(badEvent,QtCore.Qt.MatchExactly)
				if len(bEIdx) > 0:
					for B in bEIdx:
						row = self.listHigher.row(B)
						self.listHigher.item(row).setBackground(QtCore.Qt.magenta)
				bEIdx = self.listLower.findItems(badEvent,QtCore.Qt.MatchExactly)
				if len(bEIdx) > 0:
					for B in bEIdx:
						row = self.listLower.row(B)
						self.listLower.item(row).setBackground(QtCore.Qt.magenta)
				
				#Search for bad relationships in the strat relationship dropdown.
				for otherBadEvent in relListError:
					#We don't know which is older or younger, so just make two and search.
					badstring1 = str(badEvent+', '+otherBadEvent)
					badstring2 = str(otherBadEvent+', '+badEvent)
					bEIdx = self.cbxRelList.findText(badstring1)
					if bEIdx > 0:
						self.cbxRelList.model().item(bEIdx).setBackground(QtCore.Qt.magenta)
					else:
						bEIdx = self.cbxRelList.findText(badstring2)
						if bEIdx > 0:
							self.cbxRelList.model().item(bEIdx).setBackground(QtCore.Qt.magenta)
					
		else:
			#No Errors! Repaint white
			self.getRelsForCurEvent()
			if len(self.stratErrorList) > 0:
				sys.stderr.write("\nStratigraphic Error Resolved!\n")
				self.lblWarn.setText('')
				self.stratErrorList = []
			
	
	def updateListCts(self):
		if self.listHigher.count():
			self.lblHigherCt.setText(str(self.listHigher.count()))
		else:
			self.lblHigherCt.setText('')
			
		if self.listLower.count():
			self.lblLowerCt.setText(str(self.listLower.count()))
		else:
			self.lblLowerCt.setText('')
			
		if self.listUnassigned.count():
			self.lblUnassignedCt.setText(str(self.listUnassigned.count()))
		else:
			self.lblUnassignedCt.setText('')
		
		self.lblUnassignedCt.repaint()
		self.lblHigherCt.repaint()
		self.lblLowerCt.repaint()
		
		
	def selectCurEvent(self):
		self.lblWarn.setText('')
		#print(self.oldindex,'->',self.cbxCurEvent.currentIndex())
		#Add the old event back to the unassigned list if it's not "choose an event"
		#if self.oldindex > 0:
			#oldEventID = self.cbxCurEvent.itemText(self.oldindex)
			#self.listUnassigned.addItem(oldEventID)
			#self.listUnassigned.sortItems()
		self.getRelsForCurEvent()
		
		#If event is "choose an event..." gray stuff out and move all things back	
		if self.cbxCurEvent.currentIndex() <= 0:
			self.listHigher.selectAll()
			self.moveEvent('removeHigher')
			self.listLower.selectAll()
			self.moveEvent('removeLower')
			self.listHigher.setDisabled(True)
			self.listLower.setDisabled(True)
			self.listHigher.setStyleSheet('QListWidget {background-color: #eee;}')
			self.listLower.setStyleSheet('QListWidget {background-color: #eee;}')
			self.bAddHigher.setDisabled(True)
			self.bRemoveHigher.setDisabled(True)
			self.bAddLower.setDisabled(True)
			self.bRemoveLower.setDisabled(True)
		else:
			### Enable Users to work in workspace
			self.listHigher.setDisabled(False)
			self.listLower.setDisabled(False)
			self.listHigher.setStyleSheet('QListWidget {background-color: white;}')
			self.listLower.setStyleSheet('QListWidget {background-color: white;}')
			self.bAddHigher.setDisabled(False)
			self.bRemoveHigher.setDisabled(False)
			self.bAddLower.setDisabled(False)
			self.bRemoveLower.setDisabled(False)
			
			self.checkRelList()
			
		self.oldindex = self.cbxCurEvent.currentIndex()
		# Remove current event from unassigned list
		self.updateListCts()
	
	
	def getRelsForCurEvent(self):
		curEvent = self.cbxCurEvent.currentText()
		unassignedRels = []
		higherRels = []
		lowerRels = []
		
		for e in self.eventLib.events:
			unassignedRels.append(e.id)
		
		try:
			unassignedRels.remove(curEvent)
		
			for rel in self.relList:
				#if curEvent is lower, add other to higher list
				if curEvent == rel[0]:
					higherRels.append(rel[1])
					unassignedRels.remove(rel[1])
				#if curEvent is higher, add other to lower list
				elif curEvent == rel[1]:
					lowerRels.append(rel[0])
					unassignedRels.remove(rel[0])
		except ValueError:
			#Nothing to do, this means that current event is "choose an event"
			something = 0
				
		#Update GUI Lists
		self.listUnassigned.clear()
		self.listHigher.clear()
		self.listLower.clear()
		
		self.listUnassigned.addItems(sorted(unassignedRels,key=str.casefold))
		self.listHigher.addItems(sorted(higherRels,key=str.casefold))
		self.listLower.addItems(sorted(lowerRels,key=str.casefold))
		

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