#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v2/icetray-start
#METAPROJECT icerec/V05-01-00

# This code generates histograms for the likelihood that describes the reconstructed signal neutrino energy

import numpy as np
import src/ING as ING #contains histograms
import h5py
import os
import argparse
import sys

#First get the signal spectra

parser = argparse.ArgumentParser(description='Process some parameters.')
parser.add_argument('-m', dest='massnumber', default=0, type=int, choices = range(13), metavar='massnumber', help='mass number')

args    = parser.parse_args()
massnumber    = args.massnumber


masslist     = [100, 250, 350, 500, 750, 1000, 2500, 5000, 7500, 10000, 25000, 50000, 75000]
medmasslist  = [1, 10, 100, 1000, 10000]
lifetimelist = [1, 10, 100, 1000, 10000]
channellist  = [5, 8, 11, 13]
track        = []
spectrum_mu  = []
spectrum_amu = []
N_mu         = []
N_amu        = []

for mass in masslist:

    trax      = []
    spect_mu  = []
    spect_amu = []
    Nx_mu     = []
    Nx_amu    = []

    for medmass in medmasslist:

        traxx    = []
        spec_mu  = []
        spec_amu = []
        Nxx_mu   = []
        Nxx_amu  = []

        for lifetime in lifetimelist:

            traxxx   = []
            spe_mu   = []
            spe_amu  = []
            Nxxx_mu  = []
            Nxxx_amu = []

            for channel in channellist:

                if os.path.isfile("/home/ctoennis/analyses/standard_analysis_framework/spectra_wimpsim/wa-sunsum-m" + str(mass) + "-ch" + str(channel) + "-x" + str(medmass) + "-l" + str(lifetime) + "-009080-000000-read-earth_e.dat") and\
 mass == masslist[massnumber]:

                    spe_mu.append(ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/spectra_wimpsim/wa-sunsum-m" + str(mass) + "-ch" + str(channel) + "-x" + str(medmass) + "-l" + str(lifetime) + "-009080-000000-read-earth_mu.dat"))
                    Nxxx_mu.append(spe_mu[-1].Integral(spe_mu[-1].points[0][0],spe_mu[-1].points[-1][0]))
                    spe_amu.append(ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/spectra_wimpsim/wa-sunsum-m" + str(mass) + "-ch" + str(channel) + "-x" + str(medmass) + "-l" + str(lifetime) + "-009080-000000-read-earth_amu.dat"))
                    Nxxx_amu.append(spe_amu[-1].Integral(spe_amu[-1].points[0][0],spe_amu[-1].points[-1][0]))

                    traxxx.append(ING.H1D.Empty(0.0,200000.0,2000))

                else:

                    spe_mu.append("None")
                    Nxxx_mu.append("None")
                    spe_amu.append("None")
                    Nxxx_amu.append("None")

                    traxxx.append("None")


            spec_mu.append(spe_mu)
            spec_amu.append(spe_amu)
            Nxx_mu.append(Nxxx_mu)
            Nxx_amu.append(Nxxx_amu)
            traxx.append(traxxx)

        spect_mu.append(spec_mu)
        spect_amu.append(spec_amu)
        Nx_mu.append(Nxx_mu)
        Nx_amu.append(Nxx_amu)
        trax.append(traxx)

    spectrum_mu.append(spect_mu)
    spectrum_amu.append(spect_amu)
    N_mu.append(Nx_mu)
    N_amu.append(Nx_amu)
    track.append(trax)


nevents = 50000
nfiles  = 7000

#  Get a histogram that gives the likelihood to encounter the sun at a certain elevation in the sky

hsun = ING.H1D.FromFile("/home/ctoennis/analyses/standard_analysis_framework/Signal_Ingredients/HSun.txt")

hsun.Scale(hsun.nbin*1.0/(hsun.integral))


for filename in os.listdir("/data/ana/analyses/northern_tracks/version-002-p00/"): # loop over MC files

    if "MC.npy" in filename and "IC86" in filename: #Only take MC files

        infile = np.load("/data/ana/analyses/northern_tracks/version-002-p00/"+filename,"r")


        for entry in infile:

            erec   = entry["logE"]
            nutype = entry["trueType"]
            weight = entry["orig_OW"]
            energy = entry["trueE"]
	    azt    = entry["trueAzi"]
            dect   = entry["trueZen"]
            azr    = entry["azi"]
            decr   = entry["zen"]
	    k = massnumber

	    diff = (np.dot(np.array([np.sin(azt)*np.sin(dect),np.cos(azt)*np.sin(dect),np.cos(dect)]),np.array([np.sin(azr)*np.sin(decr),np.cos(azr)*np.sin(decr),np.cos(decr)]))) # Angular separation from the Sun

            if diff > 0.998 and energy <  masslist[k]: # only look at events in a small ROI around the Sun

                for j in range(len(medmasslist)):

                    for i in range(len(lifetimelist)):

                        if masslist[k] < medmasslist[j]:

                            continue

                        for l in range(len(channellist)):

                            if not spectrum_mu[k][j][i][l] == "None" and N_mu[k][j][i][l]>0:

                                if nutype == 68 or nutype == 14:

				    track[k][j][i][l].Fill(np.power(10.0,erec),weight*spectrum_mu[k][j][i][l].Eval(energy/masslist[k]))

				else:

				    track[k][j][i][l].Fill(np.power(10.0,erec),weight*spectrum_amu[k][j][i][l].Eval(energy/masslist[k]))


k = massnumber
#for k in range(len(masslist)):

for j in range(len(medmasslist)):

    for i in range(len(lifetimelist)):

        for l in range(len(channellist)):

            if not spectrum_mu[k][j][i][l] == "None" and track[k][j][i][l].integral > 0.0:

		track[k][j][i][l].Scale(len(track[k][j][i][l].content)*1.0/track[k][j][i][l].integral)

		track[k][j][i][l].Write("/home/ctoennis/analyses/standard_analysis_framework/northern_ing/ENERGY_track_m" + str(masslist[k]) + "_med" + str(medmasslist[j]) + "_l" + str(lifetimelist[i]) + "_ch" + str(channellist[l]) +  ".txt")
