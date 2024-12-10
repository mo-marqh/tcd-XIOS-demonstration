import copy
import glob
import netCDF4
import numpy as np
import os
import subprocess
import unittest

import xios_examples.shared_testing as xshared

this_path = os.path.realpath(__file__)
this_dir = os.path.dirname(this_path)

class TestRead4D(xshared._TestCase):
    test_dir = this_dir
    transient_inputs = ['sample_smaller.nc']
    transient_outputs = []
    executable = './read.exe'

    def test_read_4d(self):
        inputfile = self.transient_inputs[0]
        infile = inputfile.replace('.nc', '.cdl')
        
        subprocess.run(['ncgen', '-k', 'nc4', '-o', inputfile,
                        infile], cwd=self.test_dir, check=True)
        self.run_mpi_xios()
