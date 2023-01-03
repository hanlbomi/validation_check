def run_step4(threshold):

   import pandas as pd
   import numpy as np
   import os
   from astropy.io import fits


   res  = pd.read_csv('../check_img.dat', header = None, delimiter=r"\s+", skiprows=3) # 
   fil     = res[1]
   fil1    = res[1][0]          # 
   fil2    = res[1][len(fil)-1] #  
   nfil    = len(set(fil))
   nfil1   = np.sum(fil == fil1)
   nfil2   = np.sum(fil == fil2)
   inst    = res[3]
   chan    = res[4]
   name    = res[8]
   newname = res[9]
   rdn     = res[10]
   gain    = res[11]
   zero    = res[12]
   sfx     = res[13]



   for i in range(len(fil)):

      rad = [0.05, 0.10, 0.15]
      if inst[i] == 'ACS' and chan[i] == 'WFC':
         fwhm = 1.001 ; psl   = 0.025
      if inst[i] == 'WFC3' and chan[i] == 'UVIS':
         fwhm = 1.000 ; psl   = 0.05
      if inst[i] == 'ACS' and chan[i] == 'HRC':
         fwhm = 1.001  ; psl  = 0.07
      if inst[i] == 'WFC3' and chan[i] == 'IR':
         fwhm = 1.002  ; psl  = 0.2



      fwhm   = fwhm/psl
      fitrad = fwhm*1.5            
      psfrad = int(round(fwhm*6.)) 


      from mk_opt import mk_daophot
      mk_daophot()
      from mk_opt import mk_photo
      mk_photo()

      from run_daophot import run_daophot
      from run_daophot import cal_apphot

      name_crm1 = name[i]+'_crm_sci1'
      name_crm2 = name[i]+'_crm_sci2'
      name_crc1 = name[i]+'_crc_sci1'
      name_crc2 = name[i]+'_crc_sci2'

      print(name_crm1)

      run_daophot(name_crm1)
      os.system('rm -f daophot.log')
      cal_apphot(name_crm1)
      run_daophot(name_crc1)
      os.system('rm -f daophot.log')
      cal_apphot(name_crc1)

      if chan[i] == 'WFC' or chan[i] == 'UVIS':
         run_daophot(name_crm2)
         os.system('rm -f daophot.log')
         cal_apphot(name_crm2)
         run_daophot(name_crc2)
         os.system('rm -f daophot.log')
         cal_apphot(name_crc2)
