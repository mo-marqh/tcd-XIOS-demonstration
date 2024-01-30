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

        """
        for t_in in self.transient_inputs:
            os.remove('{}/{}'.format(self.test_dir, t_in))
        for t_out in self.transient_outputs:
            os.remove('{}/{}'.format(self.test_dir, t_out))

    @classmethod
    def tearDownClass(cls):
        """
        Finally, clean the build the fortran code only once for this class.

        """
        subprocess.run(['make', 'clean'], cwd=cls.test_dir)

