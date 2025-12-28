Statistics of Glycemia
========================

Compute All Indices
---------------------------

.. autofunction:: gly_toolbox_dev.calc_glu

.. hint:: The folder **DATA_DIR** is set to *./Data/* by default. The patient 'XX' sensor data should thus be placed in the folder *./Data/XX/*.

Example
^^^^^^^
Run with default parameters as

>>> calc_glu(patient='XX')

Provide a prescribed interval

>>> calc_glu(patient='XX',intervals=[["2025-01-12 08:00","2025-01-14 08:00"]])


List of indices implemented
---------------------------

MAGE
^^^^

A description and ref

.. autofunction:: gly_toolbox_dev.calc_MAGE

CONGAn
^^^^^^

A description and ref

.. autofunction:: gly_toolbox_dev.CONGAn

Other
^^^^^^

A description and ref


