#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v2/icetray-start
#METAPROJECT icerec/V05-01-00

import numpy as np
import src/ING as ING
import h5py
import os

hsun = ING.H1D.FromFile("/home/ctoennis/analyses/standard_analysis_framework/Signal_Ingredients/HSun.txt")

hsun = ING.H1D(-1.0,0.0,hsun.content[:hsun.GetBin(0.0)])

hsun.Scale(hsun.nbin*1.0/(hsun.integral))

track = ING.H1D.Empty(0.998,1.0,400)

#background is generated from time scrambled data by only usng event declination

for filename in os.listdir("/data/ana/analyses/northern_tracks/version-002-p00/"):

    if "exp.npy" in filename and "IC86" in filename: #Only use data files and IC86

	infile = np.load("/data/ana/analyses/northern_tracks/version-002-p00/"+filename,"r")

	for decr,azr,erec in zip(infile["zenith_MPEFit"],infile["azimuth_MPEFit"],infile["logE"]):


	    if np.power(10.0,erec) < 200000.0:

		for azt in map(lambda x: x*np.pi/20.0,range(1,40)):

                    for dect in [np.arccos(hsun.GetRandom()) for x in [0]*30]:

		        diff = (np.dot(np.array([np.sin(azt)*np.sin(dect),np.cos(azt)*np.sin(dect),np.cos(dect)]),np.array([np.sin(azr)*np.sin(decr),np.cos(azr)*np.sin(decr),np.cos(decr)])))

			if diff > 0.998:

			    track.Fill(diff,1.0/(40.0*30.0))
			


track.Write("/home/ctoennis/analyses/standard_analysis_framework/northern_ing/Background_WIMP.txt")
track.Scale(len(track.content)/track.integral)
track.Write("/home/ctoennis/analyses/standard_analysis_framework/northern_ing/Background_S_WIMP.txt")
