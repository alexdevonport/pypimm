"""
Contains utilities to process raw time-domain FMR signals, and to fit
them to a model for their behavior in time and as a function of applied
bias field. The estimated parameters have a physical significance
concerning the dynamic magnetic properties of the material that
created them.

.. automodule:: fitutils.bandlimit
    :members:

.. automodule:: fitutils.fit_data
    :members:

.. automodule:: fitutils.characterize_fits
    :members:

.. automodule:: fitutils.estimate_damping
    :members:

.. automodule:: fitutils.estimate_frequency
    :members:

.. automodule:: fitutils.fit_data
    :members:

.. automodule:: fitutils.fmin_simplex
    :members:

.. automodule:: fitutils.preprocess
    :members:

.. automodule:: fitutils.robust_fit
    :members:

.. automodule:: fitutils.shotgun_lsq
    :members:

.. automodule:: fitutils.error_analysis
    :members:

"""

from .bandlimit import bandlimit
from .fit_data import fit_data
from .characterize_fits import characterize_fits
from .estimate_damping import estimate_damping
from .estimate_frequency import estimate_frequency
from .fmin_simplex import fmin_simplex
from .preprocess import preprocess
from .robust_fit import robust_fit
from .shotgun_lsq import shotgun_lsq
from .error_analysis import error_analysis