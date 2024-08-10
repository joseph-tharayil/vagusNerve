# -*- coding: utf-8 -*-
import s4l_v1 as s4l
import numpy as np

import MeshModeling as mm
import Meshing as ms
import XCoreModeling as xcm
import XCore
import XCoreModeling
from s4l_v1 import Vec3

class ThisMesh():

    def __init__(self,mesh):
        self.Mesh=mesh
        self.Domains=[]
        self.Patches=[]
        
    def ExtractUnstructuredMeshPatchSurface(self,patch):
        xcm.ExtractUnstructuredMeshPatchSurface(patch)
        
        return
    
    def ExtractUnstructuredMeshSurface(self):
    
        xcm.ExtractUnstructuredMeshSurface(self.Mesh)
        
        return

    def ExtractUnstructuredMeshDomains(self,domainslist):

        xcm.ExtractUnstructuredMeshDomains(domainslist)
        
        return
        
    def CreateUnstructuredPatch(self,inputparams):

        xcm.CreateUnstructuredMeshPatch(inputparames)
    
        return

    def AppendUnstructuredMeshes(self,listmeshes):

        mesh=xcm.AppendUnstructuredMeshes(self,0.0,True)
        mesh.Name=newmeshname
        
        return mesh

class MaterialSettings():

    def __init__(self,tissuelistnames,tissuelist,sigmas):
        self.tissues=[]
        for i in xrange(len(tissuelist)):
            dic={'name':tissuelistnames[i],'entity':tissuelist[i],'sigma':sigmas[i]}
            self.tissues.append(dic)
        
        return

class MeshSliceSettings():

    def    __init__(self,adapt,gapele,targetl,mgr,minlength):
        self.Adaptive=adapt
        self.GapElements=gapele
        self.TargetLength=np.double(targetl)
        self.MaximumGrowRate=mgr
        self.MinLength=minlength
    
        return


def ExtrudeMesh(entitieslist,spline,npoints,meshname,invertmesh=0):
    
    tetsonly=1
    createspatches=1
    discretization=1 #1 uniform
    
    mesh=xcm.ExtrudeMesh(entitieslist,spline,tetsonly,invertmesh,createspatches,discretization,npoints)
    
    return
    
def CreateFacetedMesh(entities,options):

    options=xcm.GlobalUnstructuredMeshingOptions
    xcm.GenerateUnstructuredMesh(entities,options)

    return
    
        
ents=s4l.model.AllEntities()

label = 'Meshed Slice'

ent = ents[label]

create=0
extrude=not(create)

patches=ent.Patches


cnt=0
meshes = []

# Selects which extrusions are part of which electrode/ not in an electrode. must be updated for different nerve configurations

listBackg = [0,13,15]
listNerveCloseRec = [2,12]
listIns = [3,5,7,9,11]
listRec = [4,6,8,10]
listStim = [1020,1022]
listStimContact = [1021,14]
listStimIsh = [100,101,103,104]
listStimNear = []
# listBackg = []


lineList = [0,2,3,4,5,6,7,8,9,10,11,12,13]#,14,15]

for l in lineList: # Extrudes meshes

	line = 'Lines '+str(l)
	
	invert = 0
	
	if l in listRec or l in listStimContact:
		resolution = 0.025
	elif l in listIns:
		resolution = 0.05
	elif  l in listNerveCloseRec:
		resolution = 1
	else:
		if l == 0:
			resolution = 10
		else:
			resolution=10

	spline=ents[line]

	length=xcm.MeasureLength([spline])
	npoints=int(length/resolution)
	
	isGood = 0
	
	while not isGood:
		
		ExtrudeMesh(patches,spline,npoints,'MyMesh',invert)
		
		# isGood = 1

		if cnt==0:
			mesh=ents[label+'_extrusion']
			mesh.Name='Mesh_'+str(cnt)
		else:
			mesh=ents['Mesh_'+str(cnt-1)+'_extrusion']
			mesh.Name='Mesh_'+str(cnt)
			
		q = xcm.UnstructuredMeshQuality(mesh, xcm.eQualityType.kScaledJacobian)
		
		if np.mean(q) < 0:
			mesh.Delete()
			invert = 1
		else:
			isGood = 1
		
	meshes.append(mesh)
	patches=[f for f in mesh.Patches if 'top' in f.Name]

	cnt+=1
	
for i, iList in enumerate(lineList): # Lists of domains in each segment that are to be merged

	salineList = []
	epiList = []
	connList = []
	insList = []
	eleList = []
	stimList = []
	channelList = []

	listLists = [salineList, epiList, connList, insList, eleList]#, stimList]#, channelList]
	listNames = ['Saline','Epin','Intra','Ins','ElRec','ElStim','Channel']

	for j in range(len(listLists)): # Merges domains

		mesh = meshes[i]
		domains = mesh.Domains
		salineList = []
		epiList = []
		connList = []
		insList = []
		eleList = []
		stimList = []
		channelList = []
		
		for domain in domains: # Selects which domains need to be merged, based on which segment we are in. May need updating for different configurations
		
			name = domain.Name
			
			if iList in listBackg:
				salineList.append(domain)
				
			else:
			
				if 'Fascicles' in name:
					
					connList.append(domain)
					
				elif name == 'Nerve':
					epiList.append(domain)

				elif 'Contact' in name:
					if iList in listRec and '_Rec' in name:
						eleList.append(domain)
					elif iList in listIns and '_Rec' in name:
						insList.append(domain)
					elif iList in listStimContact:
						if '_Rec' in name:
							salineList.append(domain)
						else:
							channelList.append(domain)
					else:

						salineList.append(domain)

						
				elif 'Background' in name:
					salineList.append(domain)
						
				elif name == 'Silicone':
					if iList in listIns or iList in listRec:
						insList.append(domain)
					else:
						salineList.append(domain)
				elif 'Ele_' in name:
					if iList in listStim or iList in listStimContact:
						stimList.append(domain)
					elif "_in" in name and "ins" not in name:
						epiList.append(domain)
					elif "_out" in name:
						salineList.append(domain)
					elif "ins" in name:
						if iList in listIns or iList in listRec:
							insList.append(domain)
						else:
							salineList.append(domain)
					elif '_rec' in name:
						if iList in listRec:
							eleList.append(domain)
						elif iList in listIns:
							insList.append(domain)
						else:
							salineList.append(domain)
						
					
		listLists = [salineList, epiList, connList, insList, eleList]
		list = listLists[j]
		if len(list)>1:
			xcm.MergeUnstructuredMeshDomains(list)
		for d in mesh.Domains:
			if d.Name == 'Merged':
				d.Name = listNames[j]
				
xcm.AppendUnstructuredMeshes(meshes) # Appends meshes

mesh = ents['Mesh_0']

insList = []
for domain in mesh.Domains:
	print(domain.Name)
	if 'Ins' in domain.Name or 'Silicone' in domain.Name:
		insList.append(domain)
		
if len(insList) > 0:
	xcm.MergeUnstructuredMeshDomains(insList)
	for d in mesh.Domains:
		if d.Name == 'Merged':
			d.Name = 'Silicone'
			
nerveList = []
for domain in mesh.Domains:
	print(domain.Name)
	if 'Nerve' in domain.Name or 'Epi' in domain.Name:
		nerveList.append(domain)
		
xcm.MergeUnstructuredMeshDomains(nerveList)
for d in mesh.Domains:
	if d.Name == 'Merged':
		d.Name = 'Nerve'
		
salineList = []
for domain in mesh.Domains:
	if 'Background' in domain.Name or 'Saline' in domain.Name:
		salineList.append(domain)
	
xcm.MergeUnstructuredMeshDomains(salineList)
for d in mesh.Domains:
	if d.Name == 'Merged':
		d.Name = 'Saline'
# for domain in mesh.Domains:
	# if 'Background' in domain.Name:
		# domain.Delete()
		
for domain in mesh.Domains:
	if domain.Name == 'Intra':
		xcm.SplitUnstructuredMeshDomain(domain) #Splits fascicles

for domain in mesh.Domains:
	if domain.Name == 'Intra':
		domain.Name = 'Fascicles_0'
	elif 'Intra' in domain.Name:
		domain_num = domain.Name.split(' ')[-1]
		domain.Name = 'Fascicles_'+domain_num
		
fascicleList = []
for domain in mesh.Domains:
	if 'Fascicle' in domain.Name:
		fascicleList.append(domain)
	
xcm.MergeUnstructuredMeshDomains(fascicleList)
for d in mesh.Domains:
	if d.Name == 'Merged':
		d.Name = 'Fascicles'
