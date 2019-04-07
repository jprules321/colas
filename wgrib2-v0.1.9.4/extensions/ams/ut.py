#!/usr/bin/env python

"""
Unit tests based oy PyUnit.
"""

# Add parent directory to python search path
# ------------------------------------------
import sys
sys.path.insert(0,'..')

import os
from utudx import utUDX, run

#......................................................................

class ut(utUDX):

    def setUp(self):
        utUDX.setUp(self,['ams.udxt','ams/ams.udxt'])

    def test_getenv(self):
        self.ga('getenv USER HOME')
        self.assertEqual("USER",self.ga.rword(1,3))
        self.assertEqual("HOME",self.ga.rword(2,3))
        self.assertEqual(os.getenv("USER"),self.ga.rword(1,5))
        self.assertEqual(os.getenv("HOME"),self.ga.rword(2,5))
        try:
            self.ga('getenv DA4B928248580C274E')
        except:
            self.assertEqual("DA4B928248580C274E",self.ga.rword(1,3))
            self.assertEqual("<undef>",           self.ga.rword(1,5))

    def test_setenv(self):
        self.ga('setenv XXXX "moqueca de peixe"')
        self.assertEqual("XXXX",self.ga.rword(1,3))
        self.assertEqual("moqueca",self.ga.rword(1,5))
        self.assertEqual("de",self.ga.rword(1,6))
        self.assertEqual("peixe",self.ga.rword(1,7))
        self.ga('getenv XXXX')
        self.assertEqual("moqueca",self.ga.rword(1,5))
        self.assertEqual("de",self.ga.rword(1,6))
        self.assertEqual("peixe",self.ga.rword(1,7))
        
    def test_printenv(self):
        self.ga('setenv XXXX "moqueca de peixe"')
        self.ga('printenv $XXXX')
        self.assertEqual("moqueca",self.ga.rword(1,1))
        self.assertEqual("de",self.ga.rword(1,2))
        self.assertEqual("peixe",self.ga.rword(1,3))

    def test_pid(self):
        self.ga('printenv $$')
        self.assertEqual(os.getpid()<self.ga.rword(1,1),True)

    def test_runenv(self):
        self.ga('setenv EXPR ua;va;sqrt(ua*ua+va*va)')
        self.ga('runenv d $EXPR')
        self.ga('setenv EXPR ts')
        self.ga('@ d $EXPR')
        self.assertEqual("240",self.ga.rword(2,2))
        self.assertEqual("310",self.ga.rword(2,4))
        self.assertEqual("10", self.ga.rword(2,6))


#......................................................................

if __name__ == "__main__":
    run(ut)
