"""
Factory script for obtaining the desired interfaces from the various interface scripts.
"""

# from .mdengine import *
try:
    import mdengine
except ModuleNotFoundError:
    import atesa_v2.mdengine as mdengine
try:
    import batchsystem
except ModuleNotFoundError:
    import atesa_v2.batchsystem as batchsystem

def mdengine_factory(mdengine_toolkit):
    """
    Factory function for MDEngines.

    Parameters
    ----------
    mdengine_toolkit : str
        Name of the MDEngine to invoke.

    Returns
    -------
    None

    """
    # todo: Returns? How to document?

    mdengine_toolkits = {'amber': mdengine.AdaptAmber}

    if mdengine_toolkit not in mdengine_toolkits.keys():
        raise ValueError('unsupported MDEngine name: ' + mdengine_toolkit)

    return mdengine_toolkits[mdengine_toolkit]


def batchsystem_factory(batchsystem_toolkit):
    """
    Factory function for BatchSystems.

    Parameters
    ----------
    batchsystem_toolkit : str
        Name of the BatchSystem to invoke.

    Returns
    -------
    None

    """

    batchsystem_toolkits = {'slurm': batchsystem.AdaptSlurm, 'pbs': batchsystem.AdaptPBS, 'torque': batchsystem.AdaptPBS}   # torque and pbs are synonyms

    if batchsystem_toolkit not in batchsystem_toolkits.keys():
        raise ValueError('unsupported BatchSystem name: ' + batchsystem_toolkit)

    return batchsystem_toolkits[batchsystem_toolkit]
