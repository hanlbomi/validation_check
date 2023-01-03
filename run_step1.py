def run_step1(dir):

   import pandas as pd
   import numpy as np
   from   astroquery.mast import Observations
   from   astropy.io import fits
   from   stsci.skypac import pamutils   
   import os
   import glob
   import shutil
   from   astropy.io.fits import getdata


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
   sfx     = res[13]




   for i in range(len(fil)):
      print(name[i])

      pamutils.pam_from_file(   dir+''+name[i]+'_'+sfx[i]+'.fits', ('sci', 1), 'pam_sci1_chip2.fits')
      pam        = fits.getdata('pam_sci1_chip2.fits', header=False, ext=0)

      img_sci, head_sci = fits.getdata(dir+''+name[i]+'_aln.fits', header=True, ext=1)
      img_dq,  head_dq  = fits.getdata(dir+''+name[i]+'_aln.fits', header=True, ext=3)
      dum = (img_dq >= 2048)
      img_sci[dum] = -1000
      img_sci    = img_sci*pam
      img_sci    = np.float32(img_sci) 
      fits.writeto(''+name[i]+'_crm_sci1.fits', img_sci, head_sci, overwrite=True)

      img_crc, head_crc = fits.getdata(dir+''+name[i]+'_cln.fits', header=True, ext=1)
      img_crc    = img_crc*pam
      img_crc    = np.float32(img_crc) 
      fits.writeto(''+name[i]+'_crc_sci1.fits', img_crc, head_crc, overwrite=True)


      if chan[i] == 'WFC' or chan[i] == 'UVIS':
         pamutils.pam_from_file(dir+''+name[i]+'_'+sfx[i]+'.fits', ('sci', 2), 'pam_sci2_chip1.fits')
         pam        = fits.getdata('pam_sci2_chip1.fits', header=False, ext=0)

         img_sci, head_sci = fits.getdata(dir+''+name[i]+'_aln.fits', header=True, ext=4)
         img_dq,  head_dq  = fits.getdata(dir+''+name[i]+'_aln.fits', header=True, ext=6)
         dum = (img_dq >= 2048)
         img_sci[dum] = -1000
         img_sci    = img_sci*pam
         img_sci    = np.float32(img_sci) 
         fits.writeto(''+name[i]+'_crm_sci2.fits', img_sci, head_sci, overwrite=True)

         img_crc, head_crc = fits.getdata(dir+''+name[i]+'_cln.fits', header=True, ext=4)
         img_crc  = img_crc*pam
         img_crc  = np.float32(img_crc)
         fits.writeto(''+name[i]+'_crc_sci2.fits', img_crc, head_crc, overwrite=True)



