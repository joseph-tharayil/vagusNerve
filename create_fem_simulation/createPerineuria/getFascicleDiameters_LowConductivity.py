# -*- coding: utf-8 -*-
import s4l_v1 as s4l
import numpy as np
from s4l_v1 import Vec3
import XCoreModeling as xcm
import ModelingPostPro
from s4l_v1 import Unit
import s4l_v1.units as units

np.random.seed(2643)

s4l.ignore_deprecation_warnings()

ents=s4l.model.AllEntities()
npoints=25 #35
opt='random'
Niter=1

create=1




for nn in range(Niter):		
		
	if create:

		# Creating the simulation
		simulation = s4l.simulation.emlf.UnstructuredElectroQsOhmicSimulation()
		s4l.document.AllSimulations.Add(simulation)
		
		cnt=0
		for ent in ents['FasciclesSplit'].Entities:
		
			#slice=xcm.CreatePlanarSlice(ent,Vec3(0,0,0),Vec3(0,0,1),1)
			slice=xcm.CreatePlanarSlice(ent,Vec3(5.7,-5.5,10),Vec3(0,0,-0.363),1)
			xcm.CoverWireBody(slice)
			slice.Name='tempslice_'+ent.Name

			surf=xcm.MeasureArea([slice])
			diam=2*np.sqrt(surf/np.pi)
			
			#Adding a new ThinLayerSettings
			thin_layer_settings = s4l.simulation.emlf.ThinLayerSettings()
			components = s4l.model.AllEntities()["Mesh_0"]["Patch_"+ent.Name]
			thin_layer_settings.Name = ent.Name
			thin_layer_settings.ElectricProps.ElectricConductivity = 0.00087, Unit("S/m")
			thin_layer_settings.ThicknessIsVirtual = True
			thin_layer_settings.Thickness = 0.03*diam*1e-3, units.Meters
			simulation.Add(thin_layer_settings, components)
			
		
			slice.Delete()
	
	elif density:

		ent=s4l.model.AllEntities()['Cover Loops 1']
		a=xcm.ConvertToTriangleMesh(ent)

		meshoptions=xcm.MeshingOptions()
		meshoptions.MinEdgeLength=0.01e-3
		meshoptions.EdgeLength=0.005


		xcm.RemeshTriangleMesh(a,meshoptions)

	# point_ent.Delete()
		
		
