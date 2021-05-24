# VEAM
The *Volcanic Event Age Model* runs a Monte Carlo method to model cumulative volcanic events through time. VEAM inputs are different sources of geochronological information, that might be radiometric ages, crater retention rate modeled ages, or stratigraphic relationships. VEAM enables the modeling of recurrence rate of volcanism in a region. VEAM is currently written for Python 3.7

## CODE DESCRIPTION
**VEAM-GUI.py** - User Interface to run the Volcanic Event Age Model. See runthrough video links below.

**Strat-GUI.py** - VEAM Stratigraphy Tool to help craft a stratigraphic relationship table.

**VEAM.py** - Volcanic Event Age Model command line script developed by J Wilson. See MS Thesis at:
https://scholarcommons.usf.edu/etd/6435/

#### Python Dependencies
 - Numpy
 - Scipy
 - matplotlib
 - threading
 - queue
 - PyQt5

#### Recent Updates
 - NEW FEATURE: VEAM Stratigraphy Tool has a button to generate a list of ALL problematic stratigraphic relationships.
 - Updated VEAM-GUI.py to From Python 2.7 to Python 3.7
 - Added the VEAM Stratigraphy Tool
 - New supporting error messages when a stratigraphic issue is found in VEAM-GUI.py
 - Fixed whitespace issue with loading events
 - Fixed bug for VEAM loading very large/complex stratigraphic databases.

## SUPPLEMENTAL INFORMATION
VEAM GUI runthrough, Part 1: https://youtu.be/ZKjCvC6EYFo

VEAM GUI runthrough, Part 2: https://youtu.be/IC8lSvNN548

VEAM Stratigraphy Tool runthrough: https://youtu.be/S13zHnA55fQ

## CITATIONS
If VEAM is helpful to you, please cite:

Richardson, J. A., Wilson, J. A., Connor, C. B., Bleacher, J. E., & Kiyosugi, K. (2017). Recurrence rate and magma effusion rate for the latest volcanism on Arsia Mons, Mars. Earth and Planetary Science Letters, 458, 170-178. https://doi.org/10.1016/j.epsl.2016.10.040

Wilson, J. A., "A New Volcanic Event Recurrence Rate Model and Code For Estimating Uncertainty in Recurrence Rate and Volume Flux Through Time With Selected Examples" (2016). Graduate Theses and Dissertations. https://scholarcommons.usf.edu/etd/6435

## AUTHORS
Jacob Richardson (jacob.a.richardson@nasa.gov)

James Wilson

Aurélie Germa

Chuck Connor

Koji Kiyosugi
