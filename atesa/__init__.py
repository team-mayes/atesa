"""
atesa_v2
Version 2 of Aimless Transition Ensemble Sampling and Analysis refactors the code to make it portable, extensible, and flexible.
"""

# Add imports here
#from atesa_v2.main import *
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
