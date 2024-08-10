# -*- coding: utf-8 -*-
import s4l_v1 as s4l
import numpy as np
from s4l_v1 import Vec3
import XCoreModeling as xcm
import ModelingPostPro
from s4l_v1 import Unit
import s4l_v1.units as units

s4l.ignore_deprecation_warnings()

ents=s4l.model.AllEntities()

mesh = ents['Mesh_0']

for domain in mesh.Domains:
	if 'Fascicles' in domain.Name or 'Cover' in domain.Name or 'Electrodes' in domain.Name or 'Return' in domain.Name:
		xcm.CreateUnstructuredMeshPatch(domain)