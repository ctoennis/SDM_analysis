#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v2/icetray-start
#METAPROJECT icerec/V05-01-00

import numpy as np
import os

basefolder = "/data/user/ctoennis/WIMPSIM_DATA/"

name1 = ["prod_e","prod_ae","prod_mu","prod_amu","prod_tau","prod_atau","sur_e","sur_ae","sur_mu","sur_amu","sur_tau","sur_atau","earth_e","earth_ae","earth_mu","earth_amu","earth_tau","earth_atau","s_sur_e","s_sur_ae","s_sur_mu","s_sur_amu","s_sur_tau","s_sur_atau","s_earth_e","s_earth_ae","s_earth_mu","s_earth_amu","s_earth_tau","s_earth_atau"]

for folder in os.listdir(basefolder):

    for file in os.listdir(basefolder + folder):

        if "sunsum" in file or "md" in file:

	    infile = open(basefolder + folder + "/" + file,"r")

	    count = -1

	    print(file[:-4])
       
	    for line in infile:

	    #print(line)

	        if line[0] == "#":

		    continue

		else:

		    count += 1

		    data = line.split(" ")

		    #print(count,len(data))

		    fullname = file[:-4]+"-read-"+name1[count]+".dat"

		    outfile = open("/home/ctoennis/analyses/standard_analysis_framework/spectra_wimpsim/"+fullname,"w")
	
		    for bit,bet in zip(data[1:],map(lambda x: x/200.0,range(1,201))):

		        outfile.write(str(bet)+"      " +str(bit)+ "\n")

		    outfile.close()

	    infile.close()
