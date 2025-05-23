# -*- coding: utf-8 -*-
import s4l_v1 as s4l
import s4l_v1.simulation.neuron as neuron
import s4l_v1.document as document
import s4l_v1.model as model
import s4l_v1.analysis as analysis
import numpy as np
import matplotlib.pyplot as plt
import os
from s4l_v1 import Vec3
import XCoreModeling as xcm
import XPostProcessor as xp
import XPostProcessorUI as xpui
import s4l_v1.units as units
from scipy.interpolate import interp1d
from scipy.signal import fftconvolve, butter, sosfilt

def getCurrents():

    # Get the neuron simulation
    sim = s4l.document.AllSimulations['Nr']
    results = sim.Results()
    sensors = results.keys()

    # Initialization
    icurr_list = []
    pos_node_list = []
    inde_list = []
    neurolist = []

    neur = sensors[4]

    # Gets the transmembrane current from line sensors
    isens = results[neur]['i_membrane (*)']
    isens.Update()

    # Gets the section names and geometries
    ents = s4l.model.AllEntities()
    #neuro = ents['Lines 1 [Rat Neuron 20.00um]']
    neuro = ents['Lines 1 [Rat Neuron 4.00um]']
    geo = sim.GetSectionGeometries(neuro)
    secs = sim.GetSectionNames(neuro)

    # Gets a list of the section start, centre, and end positions
    nodecenterlist = []
    startlist = []
    endlist = []
    for key in secs:
        nodecenterlist.append(1e-6 * geo[key].SegmentCenters[0])
        startlist.append(1e-6 * geo[key].Start)
        endlist.append(1e-6 * geo[key].End)

    # Defines the positions where the membrane current is calculated
    # in this case at the center of the section
    pos_node = []
    inde = []

    for i in range(len(nodecenterlist)):
        for j in range(isens.Data.Grid.NumberOfPoints):
            d = (1e-3 * isens.Data.Grid.GetPoint(j) - nodecenterlist[i]).Length()
            if d < 1e-8:
                pos_node.append(np.array(isens.Data.Grid.GetPoint(j)))
                inde.append(j)

    ict = []
    print(inde)
    for i in range(isens.Data.NumberOfSnapshots):
        icurr = []
        ic = np.squeeze(isens.Data.Field(i))
        for l in inde:
            icurr.append(ic[l])
        ict.append(icurr)
    icurr_list.append(np.squeeze(np.array(ict)))

    currents = np.array(icurr_list)[0]
    positions = np.array(pos_node)

    nonEmpty = np.where(np.any(currents, axis=0))[0]

    currents = currents[:, nonEmpty]
    positions = positions[nonEmpty]

    return positions, currents

def get_coeffs_pointSource(positions, electrodePos, sigma=1):

    distances = np.linalg.norm(positions - electrodePos, axis=-1)
    coeffs = 1 / (4 * np.pi * sigma * distances)

    return coeffs

def getElectrodePositions(electrodeName):

    ents = s4l.model.AllEntities()
    recordingElectrode = ents[electrodeName]
    box = xcm.GetBoundingBox([recordingElectrode])
    recordingBox = np.array([np.array(box[0]), np.array(box[1])])
    recordingPos = np.mean(recordingBox, axis=0) * 1e-3

    return recordingPos

def get_pointSource_signal(electrode, referenceElectrode):

    ePos = getElectrodePositions(electrode)
    refPos = getElectrodePositions(referenceElectrode)

    nodePos, currents = getCurrents()

    coeffsRec = get_coeffs_pointSource(nodePos, ePos)
    coeffsRef = get_coeffs_pointSource(nodePos, refPos)

    signal = np.matmul(currents, coeffsRec) - np.matmul(currents, coeffsRef)
	
    print(np.where(signal<-1e-15))

    return signal

def getFlux():

    simulation = document.AllSimulations["LF"]
    simulation_extractor = simulation.Results()

    inputs = []
    model_to_grid_filter = analysis.core.ModelToGridFilter(inputs=inputs)
    model_to_grid_filter.Name = "Sphere 3"
    model_to_grid_filter.Entity = model.AllEntities()["Sphere 3"]
    model_to_grid_filter.MaximumEdgeLength = 0.0001, units.Meters
    model_to_grid_filter.UpdateAttributes()

    em_sensor_extractor = simulation_extractor["Overall Field"]
    em_sensor_extractor.FrequencySettings.ExtractedFrequency = u"All"

    inputs = [em_sensor_extractor.Outputs["J(x,y,z,f0)"], model_to_grid_filter.Outputs["Surface"]]
    field_flux_evaluator = analysis.core.FieldFluxEvaluator(inputs=inputs)
    field_flux_evaluator.UpdateAttributes()
    field_flux_evaluator.Update()

    flux = field_flux_evaluator.GetOutput(1).GetComponent(0)
    flux = np.abs(flux[0].real)

    return flux

def getTime():

    nx = 500000
    tmin = -0.03
    tmax = 0.03
    time = np.arange(tmin, tmax, (tmax - tmin) / (nx - 1))

    return time

def Scaling(deff):

    resistivity_intracellular = 0.7
    xSectionArea = np.pi * (deff / 2)**2
    resistance_intracellular = resistivity_intracellular / xSectionArea
    current_scale_factor = 1 / resistance_intracellular

    return current_scale_factor

def getVelocities(diameter):

    #d0 = 0.8e-6
    #v0 = 50e-6/(1.9e-3-1.78e-3)

    d0 = 4e-6
    v0 = 18.65
    velocity = v0 * diameter / d0

    return velocity

def loadActionPotentialShapes():

    simulation = document.AllSimulations["Nr"]
    simulation_extractor = simulation.Results()

    sensor_extractor = simulation_extractor["Overall line sensor"]
    sensor_extractor.Update()
    sensor_extractor.UpdateAttributes()
    document.AllAlgorithms.Add(sensor_extractor)

    sensor_extractor['v'].Update()

    data = sensor_extractor['v'].Data
    volts = []
    for i in range(data.NumberOfSnapshots):
        volts.append(data.Field(i)[20])

    volts = np.array(volts)
    time = np.linspace(0, 30e-3, data.NumberOfSnapshots)[:, np.newaxis]
    ap = np.hstack((time, volts))

    return ap

def FitAPShape(ap, tphi):

    tv = ap[:, 0]
    v = ap[:, 1]
    peak = tv[np.argmax(v)]
    tv -= peak

    apShapeEmpirical = v
    func = interp1d(tv, apShapeEmpirical, bounds_error=False,
                    fill_value=(apShapeEmpirical[0], apShapeEmpirical[-1]))
    Vs = func(tphi)

    # V = sosfilt(sos,Vs)
    # V[:10] = V[10]

    return Vs

def editPhiShape(phi, distance, cutoff=1e-4):

    xPos = phi[:, 0]
    xvals = xPos + distance - xPos[np.argmax(phi[:, 1])]
    phiShapeEmpirical = phi[:, 1]

    return xvals, phiShapeEmpirical

def FitPhiShape(distance, cutoff=1e-4):

    phi = ExtractPhiShape()
    np.save(r'D:\vagusNerve\SFAPValidation\phiShape.npy',phi)
    xvals = phi[:, 0]-distance
    phiShapeEmpirical = phi[:, 1]

    return interp1d(xvals, phiShapeEmpirical, bounds_error=False,
                    fill_value=(phiShapeEmpirical[0], phiShapeEmpirical[-1]))

def PhiShape(velocity, t, function):

    x = t * velocity
    out = function(x)

    return out

def ExtractPhiShape():

    simulation = document.AllSimulations["LF"]
    simulation_extractor = simulation.Results()

    inputs = []
    model_to_grid_filter = analysis.core.ModelToGridFilter(inputs=inputs)
    model_to_grid_filter.Name = "Lines 1"
    model_to_grid_filter.Entity = model.AllEntities()["Lines 1"]
    model_to_grid_filter.MaximumEdgeLength = 0.0001, units.Meters
    model_to_grid_filter.UpdateAttributes()
    document.AllAlgorithms.Add(model_to_grid_filter)

    em_sensor_extractor = simulation_extractor["Overall Field"]
    em_sensor_extractor.FrequencySettings.ExtractedFrequency = u"All"
    document.AllAlgorithms.Add(em_sensor_extractor)

    inputs = [em_sensor_extractor.Outputs["EM Potential(x,y,z,f0)"], model_to_grid_filter.Outputs["Line"]]
    field_interpolation_filter = analysis.core.FieldInterpolationFilter(inputs=inputs)
    field_interpolation_filter.UpdateAttributes()
    document.AllAlgorithms.Add(field_interpolation_filter)

    sensor = field_interpolation_filter.Outputs[0]
    sensor.Update()

    points = np.array([sensor.Data.Grid.GetPoint(i)[1] for i in range(sensor.Data.Grid.NumberOfPoints)])
    sensor.Update()
    intpot = np.real(sensor.Data.Field(0))
	
    return np.hstack((points[:, np.newaxis], intpot))

def getAnalyticSignal():

    fiberDiameter = 4e-6
    distance = 0#11*fiberDiameter*100
    stimulusDelay = 0#0.295e-3
    currentApplied = getFlux()

    aps = loadActionPotentialShapes()
    tphi = getTime()
    trueTime = np.linspace(0, 30e-3, 1000)

    timeIndices = []
    for t in trueTime:
	    timeIndices.append(np.argmin(np.abs((tphi+stimulusDelay) - t)))

    V = FitAPShape(aps, tphi)
    np.save(r'D:\vagusNerve\SFAPValidation\ap.npy',V)
	
    v = getVelocities(fiberDiameter)

    der = np.diff(V, n=2) / ((tphi[1] - tphi[0])**2)
    phiFunc = FitPhiShape(distance)
    phi = PhiShape(v, tphi, phiFunc)
    np.save(r'D:\vagusNerve\SFAPValidation\phi.npy',phi)

    cv = fftconvolve(der * Scaling(fiberDiameter) * (tphi[1] - tphi[0]) / v, phi, mode='same')
    sfap = cv / currentApplied

    return sfap[timeIndices[:-1]]

pointSourceOutput = get_pointSource_signal('Sphere 1', 'Sphere 2')
analytic = getAnalyticSignal()

np.save(r'D:\vagusNerve\SFAPValidation\analytic.npy',analytic)
np.save(r'D:\vagusNerve\SFAPValidation\pointSource.npy',pointSourceOutput)