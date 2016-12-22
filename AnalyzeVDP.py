from __future__ import division
import math
import pandas as pd
import os.path

def dropUnnecessary(df):

	drop_cols = ['Comment','Time Stamp (sec)','Status (code)','Sample Position (degrees)',
					'Bridge 1 Resistivity (Ohm-m)','Bridge 1 Excitation (uA)',
					'Bridge 2 Resistivity (Ohm-m)','Bridge 2 Excitation (uA)',
					'Bridge 3 Resistivity (Ohm-m)','Bridge 3 Excitation (uA)',
					'Bridge 4 Resistivity (Ohm-m)','Bridge 4 Excitation (uA)',
					'Bridge 1 Std. Dev. (Ohm-m)','Bridge 2 Std. Dev. (Ohm-m)',
					'Bridge 3 Std. Dev. (Ohm-m)','Bridge 4 Std. Dev. (Ohm-m)',
					'Number of Readings','Bridge 3 Resistance (Ohms)',
					'Bridge 4 Resistance (Ohms)']
	
	df = df.drop(drop_cols,axis=1)
	
	df = df.rename(columns={'Bridge 2 Resistance (Ohms)':'Ra',
							'Bridge 1 Resistance (Ohms)':'Rb'})

	return df
	
def combinePosNeg(df, fields):
	# arbitrary tolerance for field (units of Oe)
	tol = 50
	
	rxyCols = ['Temperature','Magnetic Field','Rxy']
	
	rxydf = pd.DataFrame(columns=rxyCols)
	
	for field in fields:
		negdf = df[df['Magnetic Field (Oe)'] > (-1*field-tol)]
		negdf = negdf[negdf['Magnetic Field (Oe)'] < (-1*field+tol)]
		posdf = df[df['Magnetic Field (Oe)'] > (field-tol)]
		posdf = posdf[posdf['Magnetic Field (Oe)'] < (field+tol)]
		
		posdfavg = []
		negdfavg = []
		
		for col in posdf.columns:
			posdfavg.append(posdf[col].mean(axis=0))
			
		for col in negdf.columns:
			negdfavg.append(negdf[col].mean(axis=0))
		
		negposdf = pd.concat([negdf,posdf])
		
		newrow = []
		
		newrow.append(negposdf['Temperature (K)'].mean(axis=0))
		newrow.append(0.5*(posdf['Magnetic Field (Oe)'].mean(axis=0)-negdf['Magnetic Field (Oe)'].mean(axis=0)))
		newrow.append(0.25*(posdfavg[2]+posdfavg[3]-negdfavg[2]-negdfavg[3]))
		
		rxydf = rxydf.append(pd.DataFrame([newrow],columns=rxyCols),ignore_index=True)
	
	return rxydf
	
	
def cleanFileRxy(path, fields):
	(head, extension) = os.path.splitext(path)
	newpath = head + '-cleaned' + extension
	
	df = pd.read_csv(path)
	df = dropUnnecessary(df)
	df = combinePosNeg(df, fields)
	
	df.to_csv(newpath,index=False)

	
def cleanFileRxx(path):
	'''
	This assumes that the file at path starts at the field
	labels. If your file doesn't, please delete the header.
	'''
	
	(head, extension) = os.path.splitext(path)
	newpath =  head + '-cleaned' + extension
	
	df = pd.read_csv(path)
	df = dropUnnecessary(df)
	
	# In general the resistances should not be negative for these measurements
	df = df[df['Ra'] > 0]
	df = df[df['Rb'] > 0]
	
	# There should be a way to do this without iterating through each row
	# Maybe you should work on this for future scripts
	df['Rsheet'] = 0.0
	
	for i in range(len(df.index)):
		df['Rsheet'][i] = findSheetR(df['Ra'][i],df['Rb'][i])
		#print(df['Rsheet'][i])
	
	df.to_csv(newpath,index=False)
	
def getRxxFieldSweep(path, fields):
	
	df = pd.read_csv(path)
	
	(head, extension) = os.path.splitext(path)
	newpath = head + '-cleaned' + extension
	
	tol = 50
	
	rxxCols = ['Temperature','Magnetic Field','Rxx']
	
	rxxdf = pd.DataFrame(columns=rxxCols)
	
	for field in fields:
		negdf = df[df['Magnetic Field (Oe)'] > (-1*field-tol)]
		negdf = negdf[negdf['Magnetic Field (Oe)'] < (-1*field+tol)]
		posdf = df[df['Magnetic Field (Oe)'] > (field-tol)]
		posdf = posdf[posdf['Magnetic Field (Oe)'] < (field+tol)]
		
		posdfavg = []
		negdfavg = []
		
		for col in posdf.columns:
			posdfavg.append(posdf[col].mean(axis=0))
			
		for col in negdf.columns:
			negdfavg.append(negdf[col].mean(axis=0))
		
		negposdf = pd.concat([negdf,posdf])
		
		newrow = []
		
		newrow.append(negposdf['Temperature (K)'].mean(axis=0))
		newrow.append(0.5*(posdf['Magnetic Field (Oe)'].mean(axis=0)-negdf['Magnetic Field (Oe)'].mean(axis=0)))
		newrow.append(0.5*(posdfavg[4]+negdfavg[4]))
		
		rxxdf = rxxdf.append(pd.DataFrame([newrow],columns=rxxCols),ignore_index=True)
	
	rxxdf.to_csv(newpath,index=False)
	
def findSheetR(a, b):
	error = 0.0005
	z = 2.0 * math.log1p(2) / (math.pi * (a + b))
	
	zprev = z
	
	y = 1.0 / math.exp(math.pi * zprev * a) + 1.0 / math.exp(math.pi * zprev * b)
	
	z = zprev - ((1.0 - y) / math.pi) / (a / math.exp(math.pi * zprev * a) + b / math.exp(math.pi * zprev * b))

	
	
	while ((z - zprev)/z > error):
		zprev = z
		y = 1.0 / math.exp(math.pi * zprev * a) + 1.0 / math.exp(math.pi * zprev * b)
		z = zprev - ((1.0 - y) / math.pi) / (a / math.exp(math.pi * zprev * a) + b / math.exp(math.pi * zprev * b))
	
	#print(1/z)
	return float(1.0 / z)