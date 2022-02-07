#!/usr/bin/env python
# -*- coding: utf-8 -*-

from obspy.taup.taup_create import build_taup_model
from obspy.taup.tau import TauPyModel
import numpy as np
from math import pi

"""This script helps you to esplore DeltaS-P space, based on magnitude, event-station distance, 
hypocentral depth, and stress drop"""


def kilometer2degrees(kilometer, radius=6371):
    """
    Convenience function to convert kilometers to degrees assuming a perfectly
    spherical Earth.
    :type kilometer: float
    :param kilometer: Distance in kilometers
    :type radius: int, optional
    :param radius: Radius of the Earth used for the calculation.
    :rtype: float
    :return: Distance in degrees as a floating point number.
    """
    pi = 3.14159265359
    return kilometer / (2.0 * radius * pi / 360.0)


#magnitude range:
Ml0=np.arange(2.0, 2.5, 0.5)

#distance range_
distancekm=range(20,21)

#iHypocentral depth_range:
depth=np.arange(4.5, 5, 1)

#stress drop in N/m2 (3MPa = 3e+6 Pa, 1Pa = 1 N/m2)
sdrop=3e+6


for Ml in Ml0:
    for depth1 in depth:
        for distancekm1 in distancekm:
            
            #equation used for Mo-Mw conversion in this study (from Senobary and Funning 2019)
            Mw=Ml
            Mo=10**((1.6)* Mw + 15.8) #Mo in dyne cm
            
            #dyne cm to Nm /10000000
            r3=(7*Mo/10000000)/(16*sdrop)
                              
            r=pow(r3, 1/3)

            tvel =  "./ncmodel.tvel"
            build_taup_model(tvel)
            model = TauPyModel(model="ncmodel")

            deg1 = kilometer2degrees(distancekm1)

            arrivals1P = model.get_travel_times(source_depth_in_km=depth1,distance_in_degree=deg1, phase_list=["p","P"])
            arrivals1S = model.get_travel_times(source_depth_in_km=depth1,distance_in_degree=deg1, phase_list=["s","S"])

            arr1P = arrivals1P[0]
            arr1S = arrivals1S[0]

            p1=arr1P.time
            s1=arr1S.time

            distanceM=float(r)
            distancekm2=distancekm1+(distanceM/1000)
            deg2 = kilometer2degrees(distancekm2)

            arrivals2P = model.get_travel_times(source_depth_in_km=depth1,distance_in_degree=deg2, phase_list=["p","p"])
            arrivals2S = model.get_travel_times(source_depth_in_km=depth1,distance_in_degree=deg2, phase_list=["s","S"])

            arr2P = arrivals2P[0]
            arr2S = arrivals2S[0]

            p2=arr2P.time
            s2=arr2S.time

            p12=p1-p2
            s12=s1-s2

            deltasp=s12-p12
            print("magnitude, sdrop(MPa), source_radius(m), dist(m), depth(km), epi_dist(km), s-p_threshold(s): ", Ml, sdrop/1000000, r, distanceM, depth1, distancekm1, deltasp)




