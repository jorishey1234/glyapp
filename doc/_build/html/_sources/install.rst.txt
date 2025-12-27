Installation
===============

.. autofunction:: .

Without python installation
---------------------------

Go to https://jupyter.org/try-jupyter/lab/

Create a folder Data where to upload glycemic sensors data from each patient

Create a folder Patients where the results will be saved for each patient

Make a new notebook

Copy the following code
::
	import sys, os, re, pyodide, numpy as np, pandas as pd, scipy.ndimage
	from datetime import date,time,timedelta,datetime
	import matplotlib.pyplot as plt, matplotlib.dates as mdates
	url='https://raw.githubusercontent.com/jorishey1234/glyapp/refs/heads/main/gly_toolbox_dev.py'
	exec(pyodide.http.open_url(url).read())



Without python installation
---------------------------
