import copy
import glob
import netCDF4
import numpy as np
import os
import subprocess
import unittest

this_path = os.path.realpath(__file__)
this_dir = os.path.dirname(this_path)


class _TestCase(unittest.TestCase):
    """
    UnitTest class to contain tests,
    1 test case function per input `.cdl` file

    """
    test_dir = this_dir
    transient_inputs = []
    transient_outputs = []
    rtol = 5e-03

    @classmethod
    def setUpClass(cls):
        """
        First, build the fortran code only once for this class.

        """
        subprocess.run(['make', 'clean'], cwd=cls.test_dir)
        subprocess.run(['make'], cwd=cls.test_dir)

    def tearDown(self):
        """
        After each test function, remove the input and output netCDF files.

        Use environment variable 'files' to avoid clean up of transient files
        note; this can cause multi-test classes to fail with ncgen errors, use
        for single test functions only.
        """
        for t_in in self.transient_inputs:
            rf = '{}/{}'.format(self.test_dir, t_in)
            if os.path.exists(rf) and not os.environ.get("files"):
                os.remove(rf)
        for t_out in self.transient_outputs:
            rf = '{}/{}'.format(self.test_dir, t_out)
            if os.path.exists(rf) and not os.environ.get("files"):
                os.remove(rf)

    @classmethod
    def tearDownClass(cls):
        """
        Finally, clean the build the fortran code only once for this class.

        Use environment variable 'logs' to avoid clean up, e.g. to keep logs
        """
        if not os.environ.get('logs'):
            subprocess.run(['make', 'clean'], cwd=cls.test_dir)


    @classmethod
    def make_a_resample_test(cls, inf):
        """
        this function makes a test case and returns it as a test function,
        suitable to be dynamically added to a TestCase for running.

        """
        # always copy for value, don't pass by reference.
        infile = copy.copy(inf)
        # expected by the fortran XIOS resample program's main.xml
        inputfile = cls.transient_inputs[0]
        outputfile = cls.transient_outputs[0]
        def test_resample(self):
            # create a netCDF file from the `.cdl` input
            subprocess.run(['ncgen', '-k', 'nc4', '-o', inputfile,
                            infile], cwd=cls.test_dir)
            # run the compiled Fortran XIOS programme
            run_cmd = ['mpiexec', '-n', '3', './resample.exe', ':',
                            '-n', '1', './xios_server.exe']
            if os.environ.get('MPI_FLAVOUR', '') == 'openmpi':
                # use hwthread for github ubuntu runner
                # but only for known openMPI (set by env var)
                run_cmd = run_cmd[0:1] + ['--use-hwthread-cpus'] + run_cmd[1:]
            print(' '.join(run_cmd))

            subprocess.run(run_cmd, cwd=cls.test_dir)
            # load the result netCDF file
            runfile = '{}/{}'.format(cls.test_dir, outputfile)
            assert(os.path.exists(runfile))
            rootgrp = netCDF4.Dataset(runfile, 'r')
            # read data from the resampled, expected & diff variables
            diff = rootgrp['resampled_minus_resample'][:]
            expected = rootgrp['resample_data'][:]
            result = rootgrp['resampled_data'][:]
            # prepare message for failure
            msg = ('the expected resample data array\n {exp}\n '
                   'differs from the resampled data array\n {res} \n '
                   'with diff \n {diff}\n')
            msg = msg.format(exp=expected, res=result, diff=diff)
            if np.any(diff):
                # print message for fail case,
                # as expected failures do not report msg.
                print(msg)
            # assert that all of the `diff` varaible values are zero
            # self.assertTrue(not np.any(diff), msg=msg)
            self.assertTrue(np.allclose(result, expected, rtol=cls.rtol), msg=msg)
        return test_resample
