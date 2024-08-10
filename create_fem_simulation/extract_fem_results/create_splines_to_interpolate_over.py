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

meanFiberNum = 1

Nerve = ents['Nerve']
nerveBox = xcm.GetBoundingBox([Nerve])
nerveCenter = np.mean([nerveBox[0][0],nerveBox[1][0]])

for nn in range(Niter):		
		
	if create:

		axon_ent=s4l.model.CreateGroup('Axons_Folder_'+str(nn))
		# point_ent=s4l.model.CreateGroup('Points_Folder'+str(nn))

		# Creating the simulation
		# simulation = s4l.simulation.emlf.UnstructuredElectroQsOhmicSimulation()
		# s4l.document.AllSimulations.Add(simulation)
		
		cnt=0
		for j, ent in enumerate( ents['FasciclesSplit'].Entities):
		
			slice=xcm.CreatePlanarSlice(ent,Vec3(6,-5.1,0),Vec3(0,0,1),1)
			#slice=xcm.CreatePlanarSlice(ent,Vec3(0,0,1),Vec3(0,0,-0.363),1) human
			xcm.CoverWireBody(slice)
			slice.Name='tempslice_'+ent.Name

			surf=xcm.MeasureArea([slice])
			diam=2*np.sqrt(surf/np.pi)
			print ('Surface: ', surf, diam)
			
			# Adding a new ThinLayerSettings
			# thin_layer_settings = s4l.simulation.emlf.ThinLayerSettings()
			# components = s4l.model.AllEntities()["Mesh_ 0"]["Patch_"+ent.Name.split(' ')[0]+ent.Name.split(' ')[-1]]
			# thin_layer_settings.Name = ent.Name
			# thin_layer_settings.ElectricProps.ElectricConductivity = 0.00087, Unit("S/m")
			# thin_layer_settings.ThicknessIsVirtual = False
			# thin_layer_settings.Thickness = 0.03*diam*1e-3, units.Meters
			# simulation.Add(thin_layer_settings, components)
			
			entity_bounds = xcm.GetBoundingBox([ent])
			fascCenter = np.mean([entity_bounds[0][0],entity_bounds[1][0]])
			
			nodal_distance = 1.2 #mm
			L = (entity_bounds[1][2]-entity_bounds[0][2]) - nodal_distance
			zz_min = entity_bounds[0][2] + nodal_distance/100.0 # Minimum Offset Values, also adds another small offset in case there are any boundary condition problems
			zz_max = entity_bounds[0][2] + nodal_distance # Maximum Offset Values

			# entfold=s4l.model.CreateGroup('Points_'+ent.Name)

			
			faceter = s4l.analysis.core.ModelToGridFilter()
			faceter.Entity = slice
			
			if opt=='random':
			
				numFibers = int(np.random.normal(meanFiberNum,0))

				faceter.MaximumEdgeLength = 0.05*diam
				faceter.Update(0)
				target_grid = faceter.GetOutput(0)
				#print target_grid.NumberOfPoints,faceter.MaximumEdgeLength

				r=np.uint64(np.floor(target_grid.NumberOfCells*np.random.rand(numFibers)))
				
				for i in range(numFibers):
				
					p=1e3*target_grid.GetCellCenter(r[i]); x=p[0];y=p[1]
				
					z = zz_min # np.random.uniform(zz_min,zz_max)
					

					
					a = s4l.model.CreateSpline([Vec3(x,y,z),Vec3(x,y,z+L)])
					a.Name = str(j) + '_Spline_' + str(i)+'_Random'+str(nn)
				


					axon_ent.Add(a)
		
		
			slice.Delete()
	
	elif density:

		ent=s4l.model.AllEntities()['Cover Loops 1']
		a=xcm.ConvertToTriangleMesh(ent)

		meshoptions=xcm.MeshingOptions()
		meshoptions.MinEdgeLength=0.01e-3
		meshoptions.EdgeLength=0.005


		xcm.RemeshTriangleMesh(a,meshoptions)

	# point_ent.Delete()
		
		
