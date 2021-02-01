"""
atesa
Python program for automating the "Aimless Transition Ensemble Sampling and Analysis" (ATESA) aimless shooting workflow on PBS/TORQUE or Slurm.
"""

# Add imports here
#from atesa.main import *
from . import main
from . import factory
from . import configure
from . import mdengine
from . import process
from . import interpret
from . import jobtype
from . import batchsystem
from . import utilities
from . import taskmanager
# from .configure import *
# from .factory import *
# from .mdengine import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = get_versions()['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
