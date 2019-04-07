#!/usr/bin/env python

"""
    Unit tests based on PyUnit.
"""

# Add parent directory to python search path
# ------------------------------------------
import sys
sys.path.insert(0,'..')

import os
import unittest

from utudx import utUDX, run

#......................................................................

class ut(utUDX):

    def setUp(self):
        utUDX.setUp(self,['shfilt.udxt','shfilt/shfilt.udxt'])

#--

    def test_shfilt(self):
        self._CheckCint('ps',1,4,   700,1000, 50)
        self._CheckCint('ps',5,10, -200,1000,100)
        self._CheckCint('ps',11,21,-100,1000,100)

    def test_power(self):
        self.ga('set gxout stat')
        self.ga('d sh_power(ts)')
        self.assertEqual(-1,int(self.ga.rword(2,4)))
        self.assertEqual(-999,int(self.ga.rword(4,4)))
        self.assertEqual(46,int(self.ga.rword(7,8)))
        self.assertEqual(-20000,int(self.ga.rword(9,5)))
        self.assertEqual(200000,int(self.ga.rword(9,6)))
        self.assertEqual(20000,int(self.ga.rword(9,7)))

#--
#            Useful Internal Methods for Writing Tests

    def _CheckCint(self,name,n1,n2,cmin,cmax,cint):
        """
        Check contour intervals during display.
        """
        self.ga('clear')
        self.ga('display sh_filt(%s,%d,%d)'%(name,n1,n2))
        self.assertEqual(cmin,int(self.ga.rword(1,2)))
        self.assertEqual(cmax,int(self.ga.rword(1,4)))
        self.assertEqual(cint,int(self.ga.rword(1,6)))

#......................................................................

if __name__ == "__main__":
    run(ut)

