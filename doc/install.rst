Installation
============

Without installing python
---------------------------

Go to https://jupyter.org/try-jupyter/lab/

Make a new notebook, and copy the following code in the first section. Then run the section with Shift+Enter keys 
::
	import sys, os, re, pyodide, numpy as np, pandas as pd, scipy.ndimage
	from datetime import date,time,timedelta,datetime
	import matplotlib.pyplot as plt, matplotlib.dates as mdates
	url='https://raw.githubusercontent.com/jorishey1234/glyapp/refs/heads/main/gly_toolbox_dev.py'
	exec(pyodide.http.open_url(url).read())


Start a new section, clicking on +, then if running for the first time, prepare the local environment folder structure with the command :

>>> init_environment()

Then you should be able to run glycemic statistics

>>> calc_glu('XX')

or get a plot on a week interval

>>> plot_patient('XX')

Any other patient data can now be added by drag/drop a PATIENT_NAME folder in the DATA_DIR folder in try-jupyter.


With python (v3.12) installed locally
-------------------------------------

Make sure you have the following packages installed. For instance with pip
::
	pip install sys, os, re, pyodide, numpy, pandas, scipy, datetime, matplotlib

then, in python run 

>>> from glyapp import *
>>> init_environment()
>>> calc_glu('XX')
>>> plot_patient('XX')


GlyApp File Structure
----------------------

Default Environment Variables
^^^^^^^^^^^^^^^^^^^^^^
   
.. autodata:: gly_toolbox_dev.DATA_DIR
.. autodata:: gly_toolbox_dev.RESULT_DIR

Setting Environment Folders and Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: gly_toolbox_dev.init_environment

.. autofunction:: gly_toolbox_dev.make_synthetic_gly_data


