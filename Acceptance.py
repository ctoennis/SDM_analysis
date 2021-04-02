#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v2/icetray-start
#METAPROJECT icerec/V05-01-00
import numpy as np
import src/ING as ING
import h5py
import os

masslist     = [100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000, 25000, 50000, 75000]
medmasslist  = [1, 10, 100, 1000, 10000]
lifetimelist = [1, 10, 100, 1000, 10000]
channellist  = [5, 8, 11, 13]
track        = []
spectrum_e   = []
spectrum_ae  = []
spectrum_mu  = []
spectrum_amu = []
N_e          = []
N_ae         = []
N_mu         = []
N_amu        = []

for mass in masslist:

    trax      = []
    spect_e   = []
    spect_ae  = []
    spect_mu  = []
    spect_amu = []
    Nx_e      = []
    Nx_ae     = []
    Nx_mu     = []
    Nx_amu    = []

    for medmass in medmasslist:

        traxx    = []
	spec_e   = []
	spec_ae  = []
	spec_mu  = []
	spec_amu = []
	Nxx_e    = []
	Nxx_ae   = []
	Nxx_mu   = []
	Nxx_amu  = []

        for lifetime in lifetimelist:

	    traxxx   = []
	    spe_e    = []
	    spe_ae   = []
	    spe_mu   = []
	    spe_amu  = []
	    Nxxx_e   = []
	    Nxxx_ae  = []
	    Nxxx_mu  = []
	    Nxxx_amu = []

	    for channel in channellist:

	        if os.path.isfile("/home/ctoennis/analyses/standard_analysis_framework/spectra_wimpsim/wa-sunsum-m" + str(mass) + "-ch" + str(channel) + "-x" + str(medmass) + "-l" + str(lifetime) + "-009080-000000-read-earth_e.dat"):

#		    print("wa-sunsum-m" + str(mass) + "-ch" + str(channel) + "-x" + str(medmass) + "-l" + str(lifetime) + "-009080-000000-read-earth_e.dat " + "exists")

	            spe_e.append(ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/spectra_wimpsim/wa-sunsum-m" + str(mass) + "-ch" + str(channel) + "-x" + str(medmass) + "-l" + str(lifetime) + "-009080-000000-read-earth_e.dat"))
		    Nxxx_e.append(spe_e[-1].Integral(spe_e[-1].points[0][0],spe_e[-1].points[-1][0]))
		    spe_ae.append(ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/spectra_wimpsim/wa-sunsum-m" + str(mass) + "-ch" + str(channel) + "-x" + str(medmass) + "-l" + str(lifetime) + "-009080-000000-read-earth_ae.dat"))
		    Nxxx_ae.append(spe_ae[-1].Integral(spe_ae[-1].points[0][0],spe_ae[-1].points[-1][0]))
		    spe_mu.append(ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/spectra_wimpsim/wa-sunsum-m" + str(mass) + "-ch" + str(channel) + "-x" + str(medmass) + "-l" + str(lifetime) + "-009080-000000-read-earth_mu.dat"))
		    Nxxx_mu.append(sum([x[1] for x in spe_mu[-1].points])/len(spe_mu[-1].points))
		    spe_amu.append(ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/spectra_wimpsim/wa-sunsum-m" + str(mass) + "-ch" + str(channel) + "-x" + str(medmass) + "-l" + str(lifetime) + "-009080-000000-read-earth_amu.dat"))

		    Nxxx_amu.append(sum([x[1] for x in spe_amu[-1].points])/len(spe_amu[-1].points))


		else:
		    
		    spe_e.append("None")
		    Nxxx_e.append("None")
		    spe_ae.append("None")
                    Nxxx_ae.append("None")
		    spe_mu.append("None")
                    Nxxx_mu.append("None")
		    spe_amu.append("None")
                    Nxxx_amu.append("None")

		traxxx.append(0)
		

	    spec_e.append(spe_e)
	    spec_ae.append(spe_ae)
	    spec_mu.append(spe_mu)
	    spec_amu.append(spe_amu)
	    Nxx_e.append(Nxxx_e)
	    Nxx_ae.append(Nxxx_ae)
	    Nxx_mu.append(Nxxx_mu)
	    Nxx_amu.append(Nxxx_amu)
	    traxx.append(traxxx)

	spect_e.append(spec_e)
	spect_ae.append(spec_ae)
	spect_mu.append(spec_mu)
	spect_amu.append(spec_amu)
	Nx_e.append(Nxx_e)
	Nx_ae.append(Nxx_ae)
	Nx_mu.append(Nxx_mu)
	Nx_amu.append(Nxx_amu)
	trax.append(traxx)

    spectrum_e.append(spect_e)
    spectrum_ae.append(spect_ae)
    spectrum_mu.append(spect_mu)
    spectrum_amu.append(spect_amu)
    N_e.append(Nx_e)
    N_ae.append(Nx_ae)
    N_mu.append(Nx_mu)
    N_amu.append(Nx_amu)
    track.append(trax)

Acc_casc    = ING.H1D.Empty(-1.0*np.pi,np.pi,180)
Acc_track   = ING.H1D.Empty(-1.0*np.pi,np.pi,180)


hsun = ING.H1D.FromFile("/home/ctoennis/analyses/standard_analysis_framework/Signal_Ingredients/HSun.txt")

hsun.Scale(hsun.nbin*1.0/(hsun.integral))

a = 0.0
b = 0.0

live    = 0.0

AEFF = ING.H1D.Empty(0.0,10000.0,100)
PE   = ING.H1D.Empty(0.0,10000.0,100)
NORM = ING.H1D.Empty(0.0,10000.0,100)

for filename in os.listdir("/data/ana/analyses/northern_tracks/version-002-p00"):

    if "MC.npy" in filename:

	if "IC59" in filename:

	    live = 348.138

	elif "IC79" in filename:

	    live = 310.000

	elif "2011" in filename:

	    live = 342.0883333

	else:

	    live = 1773.4748263

	live *= 0.5

	infile = np.load("/data/ana/analyses/northern_tracks/version-002-p00/"+filename,"r")
		
	nevents = len(infile["orig_OW"])

#	print nevents
	for entry in infile:

            mc     = entry["orig_OW"]
            mcen   = entry["trueE"]
            mczen  = entry["trueZen"]
            mctype = entry["trueType"]
            erec   = entry["energy_truncated"]

	    if mcen < 100000.0:

		weight = mc*live*86400.0/(10000.0)

		AEFF.Fill(mcen,weight)
		PE.Fill(mcen,weight*erec)

            if abs(np.rad2deg(mczen)-90.0)<23.4 and mcen < 100000.0:
			
		weight = mc*live*86400.0*hsun.GetContent(np.cos(mczen))/(10000.0) # 86400.0 is the number of seconds per year, 100000 is the number of square cm in a square m, hsun is a hisogram with the sun's visibility    

		for k in range(len(masslist)):

		    if(mcen < masslist[k]):

			for j in range(len(medmasslist)):

		            for i in range(len(lifetimelist)):

			        for l in range(len(channellist)):

			            if not spectrum_mu[k][j][i][l] == "None" and N_mu[k][j][i][l]>0:

					print(mctype)

					if mctype == 68 or mctype == 14:

					    track[k][j][i][l] += weight*spectrum_mu[k][j][i][l].Eval(mcen/masslist[k])/((N_mu[k][j][i][l]+N_amu[k][j][i][l])*masslist[k])

					else:
				    
					    track[k][j][i][l] += weight*spectrum_amu[k][j][i][l].Eval(mcen/masslist[k])/((N_mu[k][j][i][l]+N_amu[k][j][i][l])*masslist[k])

outfile = []

for medmass in medmasslist:

    outfil = []

    for lifetime in lifetimelist:

        outfi = []

        for channel in channellist:

	    outfi.append(open("/home/ctoennis/analyses/standard_analysis_framework/northern_ing/ACC_med" + str(medmass) + "_lifetime" + str(lifetime) + "_channel" + str(channel) + ".txt","w"))

	outfil.append(outfi)

    outfile.append(outfil)

for k in range(len(masslist)):

    for j in range(len(medmasslist)):

        for i in range(len(lifetimelist)):

            for l in range(len(channellist)):

	        if not spectrum_mu[k][j][i][l] == "None":

	            outfile[j][i][l].write(str(masslist[k]) + "    " + str(track[k][j][i][l]) + "\n")
	   


AEFF.Write("/home/ctoennis/analyses/standard_analysis_framework/northern_ing/AEFF.txt")

