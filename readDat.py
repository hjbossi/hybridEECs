# readDat.py: File to read the dat hyrbid files and write to a readable file from the jet finder
# Hannah Bossi, <hannah.bossi@cern.ch>, 06/07/2023

# import the necessary packages
import numpy as np
import pandas as pd
import sys


#Particle = make_dataclass("Particle", [("eventID", int), ("weight", float), ("cross", float), ("x", float), ("y", float)])
# load the dat file
particleList = []
with open('/home/hjb34/project/hybridModel/HYBRID_Partons_NoElastic_5020_gamma_pt40_05.dat') as f:
	eventIDBegin = int(sys.argv[1])
	eventIDEnd = int(sys.argv[2])
	weight = 0
	xsec = 0
	x = 0
	y = 0
	eventID =0
	
	
	for line in f:
		px = 0
		py = 0
		pz = 0
		mass = 0
		pdg_id = 0
		label = 0
		lineArray = line.split()
		if "weight" in lineArray[0]: 
			weight = float(lineArray[1])
			xsec = lineArray[3]
			x = lineArray[5]
			y = lineArray[7]
		elif "end" in lineArray[0]: 
			eventID = eventID +1
		elif "#" in lineArray[0]: 
			continue
			#print("---- starting new event ------")
		else: 
			px = lineArray[0]
			py = lineArray[1]
			pz = lineArray[2]
			mass = lineArray[3]
			pdg_id = lineArray[4]
			label = lineArray[5]
			if eventID > eventIDBegin:
				particleList.append({'Weight': weight, 'eventID': eventID, 'x':x, 'y': y,'px' : px, 'py':py, 'pz':pz, 'mass':mass, 'pdg_id':pdg_id, 'label':label})
		if eventID > eventIDEnd: 
			break

df = pd.DataFrame(particleList,  columns = ['Weight', 'eventID', 'x', 'y', 'px', 'py', 'pz', 'mass','pdg_id', 'label'])
df.to_csv(f'HYBRID_Partons_NoElastic_5020_gamma_pt40_05_{eventIDBegin}_{eventIDEnd}.csv')
