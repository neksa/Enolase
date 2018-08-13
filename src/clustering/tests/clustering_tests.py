#!/usr/bin/env python

import unittest
from .. import common

class Test(unittest.TestCase):

    def test_one(self):
        a = 1+1
        self.assertEqual(a, 2)

