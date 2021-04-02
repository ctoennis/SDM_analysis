#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v2/icetray-start
#METAPROJECT icerec/V05-01-00

import numpy as np
import ingredients as ING
import h5py
import os

hsun = ING.H1D.FromFile("/home/ctoennis/analyses/standard_analysis_framework/Signal_Ingredients/HSun.txt")

hsun = ING.H1D(-1.0,0.0,hsun.content[:hsun.GetBin(0.0)])

hsun.Scale(hsun.nbin*1.0/(hsun.integral))

casc  = ING.H1D.Empty(-1.0,1.0,800)
track = ING.H1D.Empty(0.998,1.0,400)

#norm = [1.0/(np.sin(np.radians(x))) for x in casc.center]

for filename in os.listdir("/data/ana/analyses/northern_tracks/version-002-p00/"):

    if "exp.npy" in filename:

	infile = np.load("/data/ana/analyses/northern_tracks/version-002-p00/"+filename,"r")

	for decr,azr,erec in zip(infile["zenith_MPEFit"],infile["azimuth_MPEFit"],infile["logE"]):

#	    print(decr)

	    if np.power(10.0,erec) < 200000.0:

		for azt in map(lambda x: x*np.pi/20.0,range(1,40)):

                    for dect in [np.arccos(hsun.GetRandom()) for x in [0]*30]:


#		        print("new")
#		        print(dect)
#			print(decr)

		        diff = (np.dot(np.array([np.sin(azt)*np.sin(dect),np.cos(azt)*np.sin(dect),np.cos(dect)]),np.array([np.sin(azr)*np.sin(decr),np.cos(azr)*np.sin(decr),np.cos(decr)])))

			if diff > 0.998:
#			    print(np.arccos(diff))

			    track.Fill(diff,1.0/(40.0*30.0))
			


#casc.Write("/home/ctoennis/analyses/standard_analysis_framework/Signal_Ingredients/Background_cscd.txt")
#casc.Scale(len(casc.content)/casc.integral)
#casc.Multiply(norm)
#casc.Write("/home/ctoennis/analyses/standard_analysis_framework/Signal_Ingredients/Background_S_cscd.txt")
track.Write("/home/ctoennis/analyses/standard_analysis_framework/northern_ing/Background_WIMP.txt")
track.Scale(len(track.content)/track.integral)
#track.Multiply(norm)
track.Write("/home/ctoennis/analyses/standard_analysis_framework/northern_ing/Background_S_WIMP.txt")
