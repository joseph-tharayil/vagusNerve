import numpy as np

from scipy.interpolate import interp1d
from scipy.signal import fftconvolve, butter, sosfilt

import s4l_v1.analysis as analysis
import s4l_v1.document as document
import s4l_v1.model as model
import s4l_v1.units as units
from s4l_v1 import ReleaseVersion
from s4l_v1 import Unit

def getTime():

    nx=500000

    tmin=-.01 # In s
    tmax=.01 # In s
    time=np.arange(tmin,tmax,(tmax-tmin)/(nx-1))

    return time

def Scaling(deff): # Diameter dependent scaling of transmembrane currents

    resistivity_intracellular = 0.7  # ohm meters

    xSectionArea = np.pi * (deff/2)**2

    resistance_intracellular = resistivity_intracellular/xSectionArea

    current_scale_factor = 1/(resistance_intracellular)

    return current_scale_factor

def getVelocities(diameter):

    d0 = 20*1e-6 # Diameters of myelinated fibers used to calculate velocities

    v0 = 97.6#56
    velocity = v0 * diameter/d0

    return velocity

def loadActionPotentialShapes():

    ReleaseVersion.set_active(ReleaseVersion.version8_0)
    
    # Creating the analysis pipeline
    # Adding a new SimulationExtractor
    simulation = document.AllSimulations["Nr"]
    simulation_extractor = simulation.Results()

    # Adding a new SensorExtractor
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
    time = np.linspace(0,2e-3,data.NumberOfSnapshots)[:,np.newaxis]

    ap = np.hstack((time,volts))

    return ap

def FitAPShape(ap,tphi): # Interpolates AP shape for a given AP

    tv = ap[:,0]
    v = ap[:,1]

    ### Sets peak time to 0
    peak = tv[np.argmax(v)]
    tv -= peak


    apShapeEmpirical = v


    func = interp1d(tv,apShapeEmpirical,bounds_error=False,fill_value=(apShapeEmpirical[0],apShapeEmpirical[-1]))

    Vs = func(tphi)


    #### Applies low-pass filter with very high cutoff, to remove artifacts
    sos = butter(1, 20000, 'lp', fs=83333, output='sos')

    #V = sosfilt(sos,Vs)

    #V[:10] = V[10]

    return Vs


def editPhiShape(phi,distance,cutoff=1e-4):

    '''
    This function takes the recording exposure curve from S4L, shifts it to match the desired distance from stimulus to recording, and smooths it
    '''

    xPos = phi[:,0]

    xvals = xPos+distance -xPos[np.argmax(phi[:,1])] # Shift to match desired distance

    phiShapeEmpirical = phi[:,1]#(phi[:,1]-np.mean(phi[:,1]))

   ########


    return xvals, phiShapeEmpirical

def FitPhiShape(distance,cutoff=1e-4):

    '''
    This function creates an interpolation object for the recording exposure
    '''

    phi = ExtractPhiShape()


    xvals, phiShapeEmpirical = editPhiShape(phi,distance,cutoff)

    return interp1d(xvals,phiShapeEmpirical,bounds_error=False,fill_value=(phiShapeEmpirical[0],phiShapeEmpirical[-1]))

def PhiShape(velocity,t,function):

    '''
    This function stretches the recording exposure in time, based on fiber velocity
    '''

    x = t*velocity ### Defines interpolation points from time vector and fiber velocity

    out = function(x)

    return out

def ExtractPhiShape():
    ReleaseVersion.set_active(ReleaseVersion.version8_0)
    
    # Creating the analysis pipeline
    # Adding a new SimulationExtractor
    simulation = document.AllSimulations["LF"]
    simulation_extractor = simulation.Results()

    # Adding a new ModelToGridFilter
    inputs = []
    model_to_grid_filter = analysis.core.ModelToGridFilter(inputs=inputs)
    model_to_grid_filter.Name = "Lines 1"
    model_to_grid_filter.Entity = model.AllEntities()["Lines 1"]
    model_to_grid_filter.MaximumEdgeLength = 0.0001, units.Meters
    model_to_grid_filter.UpdateAttributes()
    document.AllAlgorithms.Add(model_to_grid_filter)

    # Adding a new EmSensorExtractor
    em_sensor_extractor = simulation_extractor["Overall Field"]
    em_sensor_extractor.FrequencySettings.ExtractedFrequency = u"All"
    document.AllAlgorithms.Add(em_sensor_extractor)

    # Adding a new FieldInterpolationFilter
    inputs = [em_sensor_extractor.Outputs["EM Potential(x,y,z,f0)"], model_to_grid_filter.Outputs["Line"]]
    field_interpolation_filter = analysis.core.FieldInterpolationFilter(inputs=inputs)
    field_interpolation_filter.UpdateAttributes()
    document.AllAlgorithms.Add(field_interpolation_filter)

    sensor=field_interpolation_filter.Outputs[0]
    sensor.Update()

    points=np.flip(np.array([sensor.Data.Grid.GetPoint(i)[2] for i in range(sensor.Data.Grid.NumberOfPoints)]))
    sensor.Update()
    intpot=np.real(sensor.Data.Field(0))

    return np.hstack((points[:,np.newaxis],intpot))


def getFlux():

    # Creating the analysis pipeline
    # Adding a new SimulationExtractor
    simulation = document.AllSimulations["LF"]
    simulation_extractor = simulation.Results()

    # Adding a new ModelToGridFilter
    inputs = []
    model_to_grid_filter = analysis.core.ModelToGridFilter(inputs=inputs)
    model_to_grid_filter.Name = "Block 2"
    model_to_grid_filter.Entity = model.AllEntities()["Block 2"]
    model_to_grid_filter.MaximumEdgeLength = 0.0001, units.Meters
    model_to_grid_filter.UpdateAttributes()

    # Adding a new EmSensorExtractor
    em_sensor_extractor = simulation_extractor["Overall Field"]
    em_sensor_extractor.FrequencySettings.ExtractedFrequency = u"All"

    # Adding a new FieldFluxEvaluator
    inputs = [em_sensor_extractor.Outputs["J(x,y,z,f0)"], model_to_grid_filter.Outputs["Surface"]]
    field_flux_evaluator = analysis.core.FieldFluxEvaluator(inputs=inputs)
    field_flux_evaluator.UpdateAttributes()
    field_flux_evaluator.Update()

    flux=field_flux_evaluator.GetOutput(1).GetComponent(0)
    flux=np.abs(flux[0].real)

    return flux

def getAnalyticSignal():

    ### Parameters of Sim4Life Simulation

    fiberDiameter = 20e-6

    distance = 0.06 + 0.001 - 0.001# - 2000e-6*21 - 1e-6*21 # Distance between AP initiation point and recording electrode

    #distance = 0.059
    #stimulusDelay = 0.5437# There is an offset of 0.01 ms that I can't account for--should be 0.5537  # Time to action potential in Sim4Life, in ms

    stimulusDelay = 0.1034
    currentApplied = getFlux() # Current applied in S4L at the recording electrode, in amps

    ##### Performs analytic calculation #########

    #### Calculates AP shape and velocity given fiber type and diameter

    aps = loadActionPotentialShapes() # Shape of the action potential in time
    tphi=getTime()
    trueTime = np.linspace(0,2e-3,10000)

    timeIndices = []
    for t in trueTime:
        timeIndices.append(np.argmin(np.abs((tphi+stimulusDelay*1e-3)-t)))

    V = FitAPShape(aps,tphi) # Interpolates action potential shapes over time vector

    v = getVelocities(fiberDiameter)


    ###### Performs calculation

    der = np.diff(V,n=2)/((tphi[1]-tphi[0])**2) # Second derivative of AP shape

    ## Calculates exposure function
    phiFunc = FitPhiShape(distance)
    phi = PhiShape(v,tphi,phiFunc)

    ## Calculates SFAP
    cv = fftconvolve(der*Scaling(fiberDiameter)* (tphi[1]-tphi[0])/v,phi,mode='same')
    sfap = cv/(currentApplied)

    return sfap[timeIndices]

def getNeuralSensing():
    # Creating the analysis pipeline
    # Adding a new SimulationExtractor
    simulation = document.AllSimulations["LF"]
    simulation_extractor = simulation.Results()

    # Adding a new SimulationExtractor
    simulation = document.AllSimulations["Nr"]
    simulation_extractor_2 = simulation.Results()

    # Adding a new EmSensorExtractor
    em_sensor_extractor = simulation_extractor["Overall Field"]
    em_sensor_extractor.FrequencySettings.ExtractedFrequency = u"All"

    # Adding a new SensorExtractor
    sensor_extractor = simulation_extractor_2["Overall line sensor"]

    # Adding a new NeuralSensingEvaluator
    inputs = [em_sensor_extractor.Outputs["EM Potential(x,y,z,f0)"], sensor_extractor.Outputs["i_membrane (*)"]]
    neural_sensing_evaluator = analysis.neuron_evaluators.NeuralSensingEvaluator(inputs=inputs)
    neural_sensing_evaluator.Flux = np.double(getFlux()), units.Amperes

    neural_sensing_evaluator.UpdateAttributes()
    neural_sensing_evaluator.Update()
    
    return neural_sensing_evaluator.GetOutput(0).GetComponent(0)

if __name__=='__main__':
    analytic = getAnalyticSignal()
    sensing = getNeuralSensing()
    sensing *= -1
    error = np.sum(np.abs(analytic - sensing))

    np.save('analytic.npy',analytic)
    np.save('sensing.npy',sensing)
    assert error<2.3e-6