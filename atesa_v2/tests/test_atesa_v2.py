"""
Unit and regression test for the atesa_v2 package.
"""

# Import package, test suite, and other packages as needed
import atesa_v2
import pytest
import sys
from ..configure import configure
import pytraj
import os

def test_atesa_v2_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "atesa_v2" in sys.modules

def test_init_threads_restart():
    """Tests successful unpickle and return of a restart.pkl file"""
    pass # todo: implement

def test_init_threads_new():
    """Tests successful initialization of new threads"""
    settings = configure('atesa_v2/data/atesa.config')
    allthreads = atesa_v2.init_threads(settings)
    assert len(allthreads) == 1
    assert allthreads[0].coordinates == 'init.rst7'
    assert allthreads[0].topology == 'topology.prmtop'

def test_thread_get_last_frame_amber():
    """Tests thread.get_last_frame method with md_engine = 'amber'"""
    settings = configure('atesa_v2/data/atesa.config')
    settings.topology = 'atesa_v2/tests/test_data/test.prmtop'
    settings.md_engine = 'amber'
    allthreads = atesa_v2.init_threads(settings)
    allthreads[0].jobids.append('01234')
    allthreads[0].traj_files.append('atesa_v2/tests/test_data/test.nc')
    compare_traj = pytraj.iterload('atesa_v2/tests/test_data/test.rst7', allthreads[0].topology)
    query_traj = pytraj.iterload(allthreads[0].get_last_frame(0, settings), allthreads[0].topology)
    assert query_traj.n_frames == compare_traj.n_frames
    assert pytraj.center_of_mass(query_traj) == pytest.approx(pytraj.center_of_mass(compare_traj), 1e-3)
    os.remove('atesa_v2/tests/test_data/test.nc_last_frame.rst7')

def test_check_commit():
    """Tests check_commit using a dummy coordinate file"""
    settings = configure('atesa_v2/data/atesa.config')
    settings.topology = 'atesa_v2/tests/test_data/test.prmtop'
    settings.commit_fwd = [[1, 2], [3, 4], [1.0, 1.7], ['gt', 'lt']]
    settings.commit_bwd = [[1, 2], [3, 4], [0.5, 2.0], ['lt', 'gt']]
    assert atesa_v2.check_commit('atesa_v2/tests/test_data/test.rst7', settings) == 'fwd'
    settings.commit_bwd = [[1, 2], [3, 4], [1.0, 1.7], ['gt', 'lt']]
    settings.commit_fwd = [[1, 2], [3, 4], [0.5, 2.0], ['lt', 'gt']]
    assert atesa_v2.check_commit('atesa_v2/tests/test_data/test.rst7', settings) == 'bwd'
    settings.commit_fwd = [[1, 2], [3, 4], [1.0, 1.5], ['gt', 'lt']]
    settings.commit_bwd = [[1, 2], [3, 4], [0.5, 2.0], ['lt', 'gt']]
    assert atesa_v2.check_commit('atesa_v2/tests/test_data/test.rst7', settings) == ''
    settings.commit_fwd = [[1, 2], [3, 4], [1.0, 1.5], ['t', 'lt']]
    with pytest.raises(ValueError):
        atesa_v2.check_commit('atesa_v2/tests/test_data/test.rst7', settings)
    settings.commit_fwd = [[1, 2], [3, 4], [1.0, 1.5], ['gt', 'lt']]
    settings.commit_bwd = [[1, 2], [3, 4], [1.0, 1.5], ['t', 'lt']]
    with pytest.raises(ValueError):
        atesa_v2.check_commit('atesa_v2/tests/test_data/test.rst7', settings)

def test_get_status_slurm():
    """Tests thread.get_status method with batch_system = 'slurm'"""
    pass # todo: implement

def test_cancel_job_slurm():
    """Tests thread.cancel_job method with batch_system = 'slurm'"""
    pass # todo: implement
