#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 11:44:03 2019

@author: sfranke
"""

import geopy.distance

coords_1 = (52.2296756, 21.0122287)
coords_2 = (52.406374, 16.9251681)

f = geopy.distance.vincenty(coords_1, coords_2).kilometers