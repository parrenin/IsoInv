# IsoInv

A python model to invert isocrone observations in deep polar ice sheets

# Files description

* AgeModel.py: Main code which computes the age of the ice along a given radar profile.
* Maps.py: Code which can run a set of radar lines and display various results (basal age, accumulation, melting, etc.) on maps.
* DrawBemap2.py: Simple python script to draw the Bedmap2 surface map  (inset of fig. 1 of Parrenin et al., TC, 2017).
* DrawR.py: Simple python script to draw the accumulation variations with time (fig. 2 of Parrenin et al., TC, 2017).
* Maps-marie.py and Maps-marie_2.py: Modified maps used in Cavitte et al. (TC, 2018).

# How to run IsoInv?

* To run a single radar line (in spider):

	`run age_model.py radar_line_directory`

* To run a set of radar lines (in spider):

	`run maps.py radar_lines_directory`

# What is the structure of the directory containing all radar lines?

* a directory per radar line
* parameters-AllRadarLines.py: file containing parameters identical for all radar lines
* parameters-maps.py: file containing parameters for the Maps.py code.

# What is the structure of a radar line directory?

* radar-data.txt: file containing the radar datafor the radar line (see below for a description).
* parameters.py (OPT): python file containing various pramaters for the radar line. This file is optional since these parameters can be defined in the the directory upstream for all radar lines in the parameters-AllRadarLines.py file.
* ages.txt (OPT): text file with a list of ages for the isochronal layers: column1=ages, column2=sigma of ages. You can also define this file in the directory upstream for all radar lines.

# What is the output of age_model.py?

AgeModel.py creates a set of text files containing numerical outputs and pdf files containing graphs.

* a.txt: Average accumulation rate along the radar line, as well as accumulation for earch layer (as in Cavitte et al., TC, 2018).
* G0.txt: Geothermal flux along the radar line.
* m.txt: Melting rate along the radar line.
* pprime.txt: p' parameter along the radar line.
* agebottom.txt: Various ages and vertical resolution at the bottom of the ice sheet along the radar line.
* ageisochrones.txt: Average ages of isochrones as calculated by the model.
* agehorizons.txt: Average ages of the non dated horizons as calculated by the model.
* twtt.txt: Two way travel time of a set of horizons dated by the model, to compare with radar data.

* Data.pdf: Radar data along the profile (isochrones and bedrock).
* Model.pdf: Modeled age of the ice along the radar profile, as well as observed isochrones.
* Model-steady.pdf: Modeled steady age of the ice along the radar profile, as well as observed isochrones.
* Model-confidence-interval.pdf: Standard deviation of the modeled age along the radar profile.
* Temperature.pdf: Modeled temperature along the radar profile.
* Thinning.pdf: Modeled thinning function along the radar profile.
* AccumulationHistory.pdf: Layer per layer accumulation along the radar profile.
* AgeMisfit.pdf: Misfit between modeled and observed ages along the isochrones.
