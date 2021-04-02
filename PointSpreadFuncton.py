#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v2/icetray-start
#METAPROJECT icerec/V05-01-00

#This code generates the IceCube PSF for the Secluded dark matter analysis for the different  spectra provided by N. Bell et.al. JCAP 1104:003,2011


import numpy as np
import ingredients as ING # This contains a simple implementation of 1D, 2D and 3D histograms and a way to store these to files
import h5py
import os
import argparse
import sys

parser = argparse.ArgumentParser(description='Process some parameters.')
parser.add_argument('-m', dest='massnumber', default=0, type=int, choices = range(13), metavar='massnumber', help='mass number')

args    = parser.parse_args()
massnumber    = args.massnumber

#print "massnumber is " + str(massnumber)

# First get the specctra from some ASCII files generated beforehand

#dNdE_e=[map(lambda x: x+"_"+str(lifetime)+".corr",["/home/ctoennis/analyses/SDM/Spectra/e_"+str(mass) for mass in [200,500,1000,10000]]) for lifetime in [1, 100, 1000, 10000]]
#dNdE_ae=[map(lambda x: x+"_"+str(lifetime)+".corr",["/home/ctoennis/analyses/SDM/Spectra/ae_"+str(mass) for mass in [200,500,1000,10000]]) for lifetime in [1, 100, 1000, 10000]]
#dNdE_mu=[map(lambda x: x+"_"+str(lifetime)+".corr",["/home/ctoennis/analyses/SDM/Spectra/mu_"+str(mass) for mass in [200,500,1000,10000]]) for lifetime in [1, 100, 1000, 10000]]
#dNdE_amu=[map(lambda x: x+"_"+str(lifetime)+".corr",["/home/ctoennis/analyses/SDM/Spectra/amu_"+str(mass) for mass in [200,500,1000,10000]]) for lifetime in [1, 100, 1000, 10000]]

#spectrum_e    = [map(lambda x: ING.Graph.FromFile(x),y) for y in dNdE_e] 
#spectrum_ae   = [map(lambda x: ING.Graph.FromFile(x),y) for y in dNdE_ae]
#spectrum_mu   = [map(lambda x: ING.Graph.FromFile(x),y) for y in dNdE_mu]
#spectrum_amu  = [map(lambda x: ING.Graph.FromFile(x),y) for y in dNdE_amu] 

# Now get the normalisation of the Spectra

#N_e           = [map(lambda x: x.Integral(x.points[0][0],x.points[-1][0]),y) for y in spectrum_e]
#N_ae          = [map(lambda x: x.Integral(x.points[0][0],x.points[-1][0]),y) for y in spectrum_ae]
#N_mu          = [map(lambda x: x.Integral(x.points[0][0],x.points[-1][0]),y) for y in spectrum_mu]
#N_amu         = [map(lambda x: x.Integral(x.points[0][0],x.points[-1][0]),y) for y in spectrum_amu]


masslist     = [100, 250, 350, 500, 750, 1000, 2500, 5000, 7500, 10000, 25000, 50000, 75000]
medmasslist  = [1, 10, 100, 1000, 10000]
lifetimelist = [1, 10, 100, 1000, 10000]
channellist  = [5, 8, 11, 13]
track        = []
#spectrum_e   = []
#spectrum_ae  = []
spectrum_mu  = []
spectrum_amu = []
#N_e          = []
#N_ae         = []
N_mu         = []
N_amu        = []

for mass in masslist:

    trax      = []
#    spect_e   = []
#    spect_ae  = []
    spect_mu  = []
    spect_amu = []
#    Nx_e      = []
#    Nx_ae     = []
    Nx_mu     = []
    Nx_amu    = []

    for medmass in medmasslist:

        traxx    = []
#        spec_e   = []
#        spec_ae  = []
        spec_mu  = []
        spec_amu = []
#        Nxx_e    = []
#        Nxx_ae   = []
        Nxx_mu   = []
        Nxx_amu  = []

        for lifetime in lifetimelist:

            traxxx   = []
#            spe_e    = []
#            spe_ae   = []
            spe_mu   = []
            spe_amu  = []
#            Nxxx_e   = []
#            Nxxx_ae  = []
            Nxxx_mu  = []
            Nxxx_amu = []

            for channel in channellist:

                if os.path.isfile("/home/ctoennis/analyses/standard_analysis_framework/spectra_wimpsim/wa-sunsum-m" + str(mass) + "-ch" + str(channel) + "-x" + str(medmass) + "-l" + str(lifetime) + "-009080-000000-read-earth_e.dat") and mass == masslist[massnumber]:

#                   print("wa-sunsum-m" + str(mass) + "-ch" + str(channel) + "-x" + str(medmass) + "-l" + str(lifetime) + "-009080-000000-read-earth_e.dat " + "exists")

 #                   spe_e.append(ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/spectra_wimpsim/wa-sunsum-m" + str(mass) + "-ch" + str(channel) + "-x" + str(medmass) + "-l" + str(lifetime) + "-009080-000000-read-earth_e.dat"))
#                    Nxxx_e.append(spe_e[-1].Integral(spe_e[-1].points[0][0],spe_e[-1].points[-1][0]))
#                    spe_ae.append(ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/spectra_wimpsim/wa-sunsum-m" + str(mass) + "-ch" + str(channel) + "-x" + str(medmass) + "-l" + str(lifetime) + "-009080-000000-read-earth_ae.dat"))
#                    Nxxx_ae.append(spe_ae[-1].Integral(spe_ae[-1].points[0][0],spe_ae[-1].points[-1][0]))
                    spe_mu.append(ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/spectra_wimpsim/wa-sunsum-m" + str(mass) + "-ch" + str(channel) + "-x" + str(medmass) + "-l" + str(lifetime) + "-009080-000000-read-earth_mu.dat"))
                    Nxxx_mu.append(spe_mu[-1].Integral(spe_mu[-1].points[0][0],spe_mu[-1].points[-1][0]))
                    spe_amu.append(ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/spectra_wimpsim/wa-sunsum-m" + str(mass) + "-ch" + str(channel) + "-x" + str(medmass) + "-l" + str(lifetime) + "-009080-000000-read-earth_amu.dat"))
                    Nxxx_amu.append(spe_amu[-1].Integral(spe_amu[-1].points[0][0],spe_amu[-1].points[-1][0]))

                    traxxx.append(ING.H1D.Empty(0.998,1.0,400))
#                    print spe_mu[-1].Integral(spe_mu[-1].points[0][0],spe_mu[-1].points[-1][0])
#		    print [mass,medmass,lifetime,channel]

		else:

#		    spe_e.append("None")
#                    Nxxx_e.append("None")
#                    spe_ae.append("None")
#                    Nxxx_ae.append("None")
                    spe_mu.append("None")
                    Nxxx_mu.append("None")
                    spe_amu.append("None")
                    Nxxx_amu.append("None")

                    traxxx.append("None")


#            spec_e.append(spe_e)
#            spec_ae.append(spe_ae)
            spec_mu.append(spe_mu)
            spec_amu.append(spe_amu)
#            Nxx_e.append(Nxxx_e)
#            Nxx_ae.append(Nxxx_ae)
            Nxx_mu.append(Nxxx_mu)
            Nxx_amu.append(Nxxx_amu)
            traxx.append(traxxx)

#        spect_e.append(spec_e)
#        spect_ae.append(spec_ae)
        spect_mu.append(spec_mu)
        spect_amu.append(spec_amu)
#        Nx_e.append(Nxx_e)
#        Nx_ae.append(Nxx_ae)
        Nx_mu.append(Nxx_mu)
        Nx_amu.append(Nxx_amu)
        trax.append(traxx)

#    spectrum_e.append(spect_e)
#    spectrum_ae.append(spect_ae)
    spectrum_mu.append(spect_mu)
    spectrum_amu.append(spect_amu)
#    N_e.append(Nx_e)
#    N_ae.append(Nx_ae)
    N_mu.append(Nx_mu)
    N_amu.append(Nx_amu)
    track.append(trax)

#print sys.getsizeof(spectrum_e)


# Get a histogram that gives the likelihood to encounter the sun at a certain elevation in the sky

hsun = ING.H1D.FromFile("/home/ctoennis/analyses/standard_analysis_framework/Signal_Ingredients/HSun.txt")

hsun.Scale(hsun.nbin*1.0/(hsun.integral))

#Now define the output histograms

#track  = [[ING.H1D.Empty(0.98,1.0,2000),ING.H1D.Empty(0.98,1.0,2000),ING.H1D.Empty(0.98,1.0,2000),ING.H1D.Empty(0.98,1.0,2000)],[ING.H1D.Empty(0.98,1.0,2000),ING.H1D.Empty(0.98,1.0,2000),ING.H1D.Empty(0.98,1.0,2000),ING.H1D.Empty(0.98,1.0,2000)],[ING.H1D.Empty(0.98,1.0,2000),ING.H1D.Empty(0.98,1.0,2000),ING.H1D.Empty(0.98,1.0,2000),ING.H1D.Empty(0.98,1.0,2000)],[ING.H1D.Empty(0.98,1.0,2000),ING.H1D.Empty(0.98,1.0,2000),ING.H1D.Empty(0.98,1.0,2000),ING.H1D.Empty(0.98,1.0,2000)]]

count = 0.0



for filename in os.listdir("/data/ana/analyses/northern_tracks/version-002-p00/"): # loop over MC files

#    print len(os.listdir("/data/ana/analyses/northern_tracks/version-002-p00/"))
 
    if "MC.npy" in filename and "IC86" in filename: #Only take MC files

#	with np.load("/data/ana/analyses/northern_tracks/version-002-p00/"+filename,"r") as infile:
#	count = 0.0
        infile = np.load("/data/ana/analyses/northern_tracks/version-002-p00/"+filename,"r")

#	for entry in infile:

#	    print entry["trueAzi"]

	#print len(infile)

#	for azt,dect,azr,decr,nutype,weight,energy in zip(infile["trueAzi"],infile["trueZen"],infile["azi"],infile["zen"],infile["trueType"],infile["orig_OW"],infile["trueE"]) : #loop over events taking the relevant information

	for entry in infile:


	    azt    = entry["trueAzi"]
	    dect   = entry["trueZen"]
	    azr    = entry["azi"]
	    decr   = entry["zen"]
	    nutype = entry["trueType"]
	    weight = entry["orig_OW"]
	    energy = entry["trueE"]

#	    print"\r"
#	    sys.stdout.write("\033[F")
#	    sys.stdout.write("\033[K")
#	    print 100*count/len(infile)

#	    count += 1

	    diff = (np.dot(np.array([np.sin(azt)*np.sin(dect),np.cos(azt)*np.sin(dect),np.cos(dect)]),np.array([np.sin(azr)*np.sin(decr),np.cos(azr)*np.sin(decr),np.cos(decr)]))) # Angular separation from the Sun

	    k = massnumber

	    if diff > 0.998 and energy <  masslist[k]: # only look at events in a small ROI around the Sun

		#count += 1.0

		#for k in range(len(masslist)):

               
		
                for j in range(len(medmasslist)):

                    for i in range(len(lifetimelist)):

		        if masslist[k] < medmasslist[j]:

			    continue

                        for l in range(len(channellist)):

                            if not spectrum_mu[k][j][i][l] == "None" and N_mu[k][j][i][l]>0:

                                if nutype == 68 or nutype == 14:


#				    print spectrum_mu[k][j][i][l].Eval(energy/masslist[k])
				    track[k][j][i][l].Fill(diff,weight*spectrum_mu[k][j][i][l].Eval(energy/masslist[k]))
					
				else: # take care of antineutrinos

#					    print "amu"
				    track[k][j][i][l].Fill(diff,weight*spectrum_amu[k][j][i][l].Eval(energy/masslist[k]))

        del infile
# Store the PSF in ASCII files now
k = massnumber
#for k in range(len(masslist)):

for j in range(len(medmasslist)):

    for i in range(len(lifetimelist)):

        for l in range(len(channellist)):

            if not spectrum_mu[k][j][i][l] == "None":

		track[k][j][i][l].Write("/home/ctoennis/analyses/standard_analysis_framework/northern_ing/SPSF_WIMPSIM_m" + str(masslist[k]) + "_med" + str(medmasslist[j]) + "_l" + str(lifetimelist[i]) + "_ch" + str(channellist[l]) + ".txt")
	        track[k][j][i][l].Scale(len(track[k][j][i][l].content)*1.0/track[k][j][i][l].integral)
		track[k][j][i][l].Write("/home/ctoennis/analyses/standard_analysis_framework/northern_ing/PSF_WIMPSIM_m" + str(masslist[k]) + "_med" + str(medmasslist[j]) + "_l" + str(lifetimelist[i]) + "_ch" + str(channellist[l]) + ".txt")

