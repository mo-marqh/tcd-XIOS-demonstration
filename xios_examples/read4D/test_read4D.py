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

    @unittest.skipIf(os.environ.get('MVER', '') == 'XIOS/trunk@2252',
                     "skipping for ")
    def test_read_4d(self):
        inputfile = self.transient_inputs[0]
        infile = inputfile.replace('.nc', '.cdl')
        
        subprocess.run(['ncgen', '-k', 'nc4', '-o', inputfile,
                        infile], cwd=self.test_dir, check=True)
        self.run_mpi_xios()

    @classmethod
    def setUpClass(cls):
        if os.environ.get('MVER', '').startswith('XIOS3/trunk'):
            with open(os.path.join(cls.test_dir, 'read.F90'), 'r') as ioin:
                iodef_in = ioin.read()
            # patch in XIOS3 domain access source code syntax
            in2 = 'xios_get_domain_attr("original_domain"'
            in3 = ('xios_get_domain_attr("specific_humidity::"')
            iodef_out = iodef_in.replace(in2, in3)
            with open(os.path.join(cls.test_dir, 'read.F90'), 'w') as ioout:
                ioout.write(iodef_out)
        super().setUpClass()
