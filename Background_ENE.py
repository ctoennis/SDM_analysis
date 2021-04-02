#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v2/icetray-start
#METAPROJECT icerec/V05-01-00

import numpy as np
import ingredients as ING
import h5py
import os

hsun = ING.H1D.FromFile("/home/ctoennis/analyses/standard_analysis_framework/Signal_Ingredients/HSun.txt")

hsun.Scale(hsun.nbin*1.0/(hsun.integral))

track = ING.H1D.Empty(0.0,200000.0,2000)

#norm = [1.0/(np.sin(np.radians(x))) for x in casc.center]

for filename in os.listdir("/data/ana/analyses/northern_tracks/version-002-p00/"):

    if "exp.npy" in filename:

	infile = np.load("/data/ana/analyses/northern_tracks/version-002-p00/"+filename,"r")

	for decr,erec in zip(infile["dec"],infile["logE"]):

	    if abs(decr) < 50.0*np.pi/180.0 :

		track.Fill(np.power(10.0,erec))

track.Write("/home/ctoennis/analyses/standard_analysis_framework/northern_ing/Background_WIMP_ene.txt")
track.Scale(len(track.content)/track.integral)
#track.Multiply(norm)
track.Write("/home/ctoennis/analyses/standard_analysis_framework/northern_ing/Background_S_WIMP_ene.txt")
