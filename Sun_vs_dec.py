import argparse
import numpy as np
from icecube import dataclasses,astro
import ingredients as ing

HSun = ing.H1D.Empty(-1.0,0.0,200)

TIME=dataclasses.I3Time(0)

for year in [2010,2011,2012,2013,2014,2015,2016]:

    print year

    for month in range(3,11):

        print month

        for day in range(1,32):

            if month==2 and day > 28:

                if not year in [2012,2016]:

                    continue

                elif day > 29:

                    continue

            if month in [4,6,9,11] and day==31:

                continue
            
            for hour in range(24):

                for minute in range(0,60):

                    #for second in range(0,60,10):
                        
                    TIME.set_utc_cal_date(year,month,day,hour,minute,0,0.0)
                    print minute
                    print TIME
                    pos=np.cos(astro.I3GetSunDirection(TIME).zenith)

                    if pos < 0.0:

#                        print pos

                        HSun.Fill(pos)

HSun.Write("/home/ctoennis/Analyses/Standard_Analysis_Framework/Signal_Ingredients/HSun.txt")
