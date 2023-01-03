# ISJang's  pipeline (ver 1.0, 29Nov2020)

#Run the following modules one by one (yes/no). 


# Copy images
step1          = 'no'
dir            = '../../data/'

# Run ap.phot for the diagram A
step2          = 'no'
threshold2     = 6


# Plot the diagram A
step3          = 'no'
FID3           = 2
med            = 1.08 + 0.06
wid            = 0.08 + 0.04
mbri           = 19.7 #+ 0.7
mfai           = 23.5 + 0.5
mhist          = 24.0 + 1
bs             = 0.025
cra            = 0.0
cdec           = 0.0
rad            = 0.0
cini           = 0.7
cfin           = 1.9
mini           = 26.5
mfin           = 18.0

# Run ap.phot for the CR cnadidates
step4          = 'no'
threshold4     = 10

# apply inst. mag calibration 
step5          = 'no'

# multi-dimem. PDFs 
step6          = 'no'
FID6           = 2

# validation check
step7          = 'yes'



#====================================================================
#=============== Nothing to modify below this line ==================
#====================================================================

import sys
sys.path.append("/python_lib/")

from      glob  import glob
from run_step1  import run_step1
from run_step2  import run_step2
from run_step3  import run_step3
from run_step4  import run_step4
from run_step5  import run_step5
from run_step6  import run_step6
from run_step7  import run_step7



if step1  == 'yes':
   run_step1(dir)

if step2  == 'yes':
   run_step2(threshold2)

if step3  == 'yes':
   run_step3(dir, FID3, med, wid, mbri, mfai, mhist, bs, cra, cdec, rad, cini, cfin, mini, mfin)

if step4  == 'yes':
   run_step4(threshold4)

if step5  == 'yes':
   run_step5(dir)

if step6  =='yes':
   run_step6(FID6)

if step7  =='yes':
   run_step7()
