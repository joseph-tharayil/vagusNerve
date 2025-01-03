# -*- coding: utf-8 -*-
# This script was auto-generated by Sim4Life version 7.2.1.11125

import numpy
import s4l_v1.analysis as analysis
import s4l_v1.document as document
import s4l_v1.model as model
import s4l_v1.units as units
from s4l_v1 import ReleaseVersion
from s4l_v1 import Unit

try:
	# Define the version to use for default values
	ReleaseVersion.set_active(ReleaseVersion.version7_2)

	# Creating the analysis pipeline
	# Adding a new SimulationExtractor
	simulation = document.AllSimulations["SuperFar"]
	simulation_extractor = simulation.Results()

	# Adding a new ModelToGridFilter
	inputs = []
	
	
	for i in range(39):
	
		inputs = []
		model_to_grid_filter = analysis.core.ModelToGridFilter(inputs=inputs)
	
		model_to_grid_filter.Name = str(i)+"_Spline_0_Random0"
		model_to_grid_filter.Entity = model.AllEntities()[str(i)+"_Spline_0_Random0"]
		model_to_grid_filter.MaximumEdgeLength = 1e-6, units.Meters
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

		# Adding a new ExcelExporter
		inputs = [field_interpolation_filter.Outputs["EM Potential(x,y,z,f0)"]]
		excel_exporter = analysis.exporters.ExcelExporter(inputs=inputs)
		excel_exporter.FileName = u"C:\\Users\\tharayil\\Desktop\\SimResults\\FasciclePhiBigFar\\"+str(i)+"_BetterConductivity.xlsx"
		excel_exporter.UpdateAttributes()
		document.AllAlgorithms.Add(excel_exporter)

except Exception as exc:
	import traceback
	traceback.print_exc()
	# Reset active version to default
	ReleaseVersion.reset()
	raise(exc)
