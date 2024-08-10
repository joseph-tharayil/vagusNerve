# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
# Take in the list of splines or points within the folder "Splines_Small_Areas"
# and "Splines_Big_Areas" respectively.
import s4l_v1.document as document
import s4l_v1.materials.database as database
import s4l_v1.model as model
import s4l_v1.simulation.neuron as neuron
import s4l_v1.units as units
from s4l_v1 import ReleaseVersion
from s4l_v1 import Unit
import pickle
from scipy.stats import rv_histogram
import s4l_v1 as s4l
import numpy as np
from s4l_v1 import Vec3
import XCoreModeling as xcm
import ModelingPostPro

def getDiameters(type,fascIdx,numFibers):
    

    
    maffvals = np.loadtxt(r'C:\Users\tharayil\Downloads\maffvals.csv',delimiter=',')
    meffvals = np.loadtxt(r'C:\Users\tharayil\Downloads\meffvalsSmooth.csv',delimiter=',')

    
    h = rv_histogram((maffvals[:-1,1],maffvals[:,0]))
    diams  = h.rvs(size=numFibers[0],random_state=fascIdx)
    
    h2 = rv_histogram((meffvals[:-1,1],meffvals[:,0]))
    diams2  = h2.rvs(size=numFibers[1],random_state=fascIdx)
	
    d = [diams,diams2]

    
    return d[type]


np.random.seed(2643)

s4l.ignore_deprecation_warnings()

ents=s4l.model.AllEntities()


meanFiberNum = 50

Nerve = ents['Nerve']
nerveBox = xcm.GetBoundingBox([Nerve])
nerveCenter = np.mean([nerveBox[0][0],nerveBox[1][0]])

		

axon_ent=s4l.model.CreateGroup('Axons_Folder_2')

for ent in ents['Fascicles'].Entities:

	slice=xcm.CreatePlanarSlice(ent,Vec3(6,-5.1,-0),Vec3(0,0,1),1)
	xcm.CoverWireBody(slice)
	slice.Name='tempslice_'+ent.Name

	surf=xcm.MeasureArea([slice])
	diam=2*np.sqrt(surf/np.pi)
	
	entity_bounds = xcm.GetBoundingBox([ent])

	nodal_distance = 1.2 #mm
	L = (entity_bounds[1][2]-entity_bounds[0][2]) - nodal_distance + 20
	zz_min = entity_bounds[0][2] + nodal_distance/100.0 # Minimum Offset Values, also adds another small offset in case there are any boundary condition problems
	zz_max = entity_bounds[0][2] + nodal_distance # Maximum Offset Values

	
	splinesUeff = s4l.model.CreateGroup('SplinesUeff_'+ ent.Name)

	splinesList = [splinesUeff]

	faceter = s4l.analysis.core.ModelToGridFilter()
	faceter.Entity = slice
	
	
	numFibers = int(np.random.normal(meanFiberNum,0))

	faceter.MaximumEdgeLength = 0.05*diam
	faceter.Update(0)
	target_grid = faceter.GetOutput(0)

	r=np.uint64(np.floor(target_grid.NumberOfCells*np.random.rand(numFibers)))
	
	ueffFrac = 1
	
	probs = np.array([ueffFrac])
	typeNames = ['ueff']
	
	numFibersByType = (probs*numFibers).astype(int)
	
	for type, p in enumerate(probs):
	
		numType = int(p*numFibers)
		
		if type == 0:
			numUeff = numType
	
		for i in range(numType):
		
			if type == 0:
				j = i
			else:
				j = i + numUeff
		
			p=1e3*target_grid.GetCellCenter(r[j]); x=p[0];y=p[1]
		
			z = np.random.uniform(zz_min,zz_max)
			
			
			a = s4l.model.CreateSpline([Vec3(x,y,z),Vec3(x,y,z+L)])
			a.Name = ent.Name +'_'+typeNames[type]+ '_Spline_' + str(j)
			
			b = s4l.model.CreatePoint(Vec3(x,y,z))
			b.Name = a.Name
		
			splinesList[type].Add(a)

	axon_ent.Add(splinesUeff)

	slice.Delete()



# Define the version to use for default values
ents=model.AllEntities()
    
# Creating the simulation

simulation = neuron.Simulation()
simulation.Name = "Titrate_Schild"
document.AllSimulations.Add( simulation )

for i, fold in enumerate(ents['Axons_Folder_2'].Entities):

	fascIdx = int(i/2)
	type = i%2
	
	diameters = 0.8

	for j, ent in enumerate(fold.Entities):
	
		r = diameters
	
		automatic_axon_neuron_settings = neuron.AutomaticAxonNeuronSettings()
	
		sundt_settings = model.SchildNeuronProperties()
		components = [ent]

		sundt_settings.AxonDiameter = r
		nrn = model.CreateAxonNeurons(components, sundt_settings)

			
		simulation.Add(automatic_axon_neuron_settings, nrn)
			
		

# Load the Model        
simulation.LoadModel()

neurons=simulation.GetNeuronEntities()

# for ent in neurons:
	# line_sensor_settings = simulation.AddLineSensor(ent)
	# line_sensor_settings.NumberOfSnapshots=400
	# line_sensor_settings.Quantity=1	
	# line_sensor_settings.RecordV = False
	# line_sensor_settings.RecordIMembrane = False
	# line_sensor_settings.RecordScaledIMembrane = True