#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v2/icetray-start
#METAPROJECT icerec/V05-01-00

import numpy as np
import ingredients as ING
from scipy.stats import poisson
import os

TS = []
masslist     = [100, 250, 350, 500, 750, 1000, 2500, 5000, 7500, 10000, 25000, 50000, 75000]
medmasslist  = [1, 10, 100, 1000, 10000]
lifetimelist = [1, 10, 100, 1000, 10000]
channellist  = [5, 8, 11, 13]

for mass in masslist:

    TS_dest = []

    for channel in channellist:

        TS_dust = []

        for life in lifetimelist:

            TS_dost = []

            for medmass in medmasslist:

	        TS_dist = []

	        if os.path.isfile("/home/ctoennis/analyses/standard_analysis_framework/new_TS/TS_m" + str(mass) + "-med" + str(medmass) + "-gl" + str(life) + "-ch" + str(channel) + "-ns0.txt"):

#		    print("found something")

		    for n_s in range(61):

		        if os.path.isfile("/home/ctoennis/analyses/standard_analysis_framework/new_TS/TS_m" + str(mass) + "-med" + str(medmass) + "-gl" + str(life) + "-ch" + str(channel) + "-ns"+ str(n_s)+".txt"):

		            TS_dist.append(ING.H1D.FromFile("/home/ctoennis/analyses/standard_analysis_framework/new_TS/TS_m" + str(mass) + "-med" + str(medmass) + "-gl" + str(life) + "-ch" + str(channel) + "-ns"+ str(n_s)+".txt"))

			else:

			    TS_dist.append("none")
		
		else:

		    for n_s in range(61):

		        TS_dist.append("none")

	        TS_dost.append(TS_dist)
	
	    TS_dust.append(TS_dost)
	    
	TS_dest.append(TS_dust)

    TS.append(TS_dest)

#print(TS[0][3])

pois = []

for mass in range(len(masslist)):

    pois_dest = []

    for channel in range(len(channellist)):

        pois_dust = []

        for life in range(len(lifetimelist)):

            pois_dost = []

            for medmass in range(len(medmasslist)):

	        pois_dist = []

		if os.path.isfile("/home/ctoennis/analyses/standard_analysis_framework/new_TS/TS_m" + str(masslist[mass]) + "-med" + str(medmasslist[medmass]) + "-gl" + str(lifetimelist[life]) + "-ch" + str(channellist[channel]) + "-ns0.txt"):

#		    print("found something")

                    for n_s in [0.1*x for x in range(1,601)]:

		        h_new = ING.H1D.Empty(TS[0][3][0][0][0].low, TS[0][3][0][0][0].high, TS[0][3][0][0][0].nbin)

			for n_2 in range(60):

			    if not TS[mass][channel][life][medmass][n_2] == "none":
			
				h_new.AddH1(TS[mass][channel][life][medmass][n_2],poisson.pmf(n_2,n_s))


			h_new.Write("/home/ctoennis/analyses/standard_analysis_framework/new_TS/POIS_m" + str(mass) + "-med" + str(medmass) + "-gl" + str(life) + "-ch" + str(channel) + "_ns" +str(10*n_s) + ".txt")

			pois_dist.append(h_new)

		else:

		    for n_s in [0.1*x for x in range(1,601)]:

		        pois_dist.append("none")

		pois_dost.append(pois_dist)

	    pois_dust.append(pois_dost)

	pois_dest.append(pois_dust)

    pois.append(pois_dest)

outfile = []

#print pois[3][3][190].content
for channel in range(len(channellist)):

    outfilee = []

    for life in range(len(lifetimelist)):

        outfileee = []

        for medmass in range(len(medmasslist)):
	
	    if not TS[mass][channel][life][medmass][0] == "none":

		outfileee.append(open("/home/ctoennis/analyses/standard_analysis_framework/new_sens/ns_limit-med" + str(medmass) + "-gl" + str(life) + "-ch" + str(channel) + ".txt",'w'))

	    else:

		outfileee.append("none")

	outfilee.append(outfileee)

    outfile.append(outfilee)

for mass in range(len(masslist)):

    for channel in range(len(channellist)):

        for life in range(len(lifetimelist)):

            for medmass in range(len(medmasslist)):

	        if not TS[mass][channel][life][medmass][0] == "none":

                    back  = TS[mass][channel][life][medmass][0].GetFCInterval(0.5)

		#print back

		    derp = [0,0]

		    limit = 60.0

		    for n_s in range(1,600):
		    		    
		        if not pois[mass][channel][life][medmass][n_s] == "none":

		            derp = pois[mass][channel][life][medmass][n_s].GetFCInterval(0.9)
			#	     print derp[0], back[1], n_s
			
			    if derp[0] > back[1]:

				limit = (n_s-1)*0.1
				
				break

		    if not  outfile[channel][life][medmass] == "none":

			outfile[channel][life][medmass].write(str(masslist[mass]) + "    " + str(limit) + "\n")

			print str(masslist[mass]) + "    " + str(life) + "    " + str(limit)
		    



