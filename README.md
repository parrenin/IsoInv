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

	`run AgeModel.py radar_line_directory`

* To run a set of radar lines (in spider):

	`run Maps.py radar_lines_directory`

# What is the structure of the directory containing all radar lines?

* a directory per radar line
* parameters-AllRadarLines.py: file containing parameters identical for all radar lines
* parameters-maps.py: file containing parameters for the Maps.py code.

# What is the structure of a radar line directory?

* radar-data.txt: file containing the radar datafor the radar line (see below for a description).
* parameters.py (OPT): python file containing various pramaters for the radar line. This file is optional since these parameters can be defined in the the directory upstream for all radar lines in the parameters-AllRadarLines.py file.
* ages.txt (OPT): text file with a list of ages for the isochronal layers: column1=ages, column2=sigma of ages. You can also define this file in the directory upstream for all radar lines.

