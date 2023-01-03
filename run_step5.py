def run_step5(dir):

   import pandas as pd
   import numpy as np
   import os
   from astropy.io import fits
   import astropy.wcs as wcs
   from astropy import units as u
   from astropy.coordinates import SkyCoord



   res  = pd.read_csv('../check_img.dat', header = None, delimiter=r"\s+", skiprows=3) # 
   fil     = res[1]
   fil1    = res[1][0]          # 
   fil2    = res[1][len(fil)-1] #  
   nfil    = len(set(fil))
   nfil1   = np.sum(fil == fil1)
   nfil2   = np.sum(fil == fil2)
   expo    = res[2]
   inst    = res[3]
   chan    = res[4]
   name    = res[8]
   newname = res[9]
   rdn     = res[10]
   gain    = res[11]
   zero    = res[12]
   sfx     = res[13]



   for i in range(len(fil)):

      nelex = 4096
      neley = 2048
      if chan[i] == 'HRC' or chan[i] == 'IR':
         nelex = 1024
         neley = 1024

 
      print(''+str(fil[i])+'_psf.reg')
      res  = pd.read_csv(''+str(fil[i])+'_psf.reg', header=None, skiprows=3, engine='python', sep='[(,]')
      xdrc = res[1]
      ydrc = res[2]

      hdulist = fits.open('../'+str(fil[i])+'_sci.fits')
      w0 = wcs.WCS(hdulist[0].header, hdulist)
      hdulist.close()
      ra, dec = w0.all_pix2world(xdrc, ydrc, 1)

      fmt = '%6d %9.3f %9.3f %9.3f %9.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n'



      hdulist = fits.open(dir+''+name[i]+'_'+sfx[i]+'.fits')
      w1 = wcs.WCS(hdulist['sci',1].header, hdulist)
      hdulist.close()

      xpsf, ypsf = w1.all_world2pix(ra, dec, 1) 

      res   = pd.read_csv(''+name[i]+'_crm_sci1.ap.mag', header=None, delimiter=r"\s+") 
      x_crm    = res[1]  ; y_crm    = res[2]
      m01_crm  = res[3]  ; m02_crm  = res[4]  ; m03_crm  = res[5]  ; m04_crm  = res[6]  
      m05_crm  = res[7]  ; m06_crm  = res[8]  ; m07_crm  = res[9]  ; m08_crm  = res[10] 
      m09_crm  = res[11] ; m10_crm  = res[12] ; m11_crm  = res[13] ; m12_crm  = res[14]

      res   = pd.read_csv(''+name[i]+'_crm_sci1.ap.err', header=None, delimiter=r"\s+") 
      e01_crm  = res[3]  ; e02_crm  = res[4]  ; e03_crm  = res[5]  ; e04_crm  = res[6]
      e05_crm  = res[7]  ; e06_crm  = res[8]  ; e07_crm  = res[9]  ; e08_crm  = res[10]
      e09_crm  = res[11] ; e10_crm  = res[12] ; e11_crm  = res[13] ; e12_crm  = res[14]

      res   = pd.read_csv(''+name[i]+'_crc_sci1.ap.mag', header=None, delimiter=r"\s+") 
      x_crc    = res[1]  ; y_crc    = res[2]
      m01_crc  = res[3]  ; m02_crc  = res[4]  ; m03_crc  = res[5]  ; m04_crc  = res[6]  
      m05_crc  = res[7]  ; m06_crc  = res[8]  ; m07_crc  = res[9]  ; m08_crc  = res[10] 
      m09_crc  = res[11] ; m10_crc  = res[12] ; m11_crc  = res[13] ; m12_crc  = res[14]

      res   = pd.read_csv(''+name[i]+'_crc_sci1.ap.err', header=None, delimiter=r"\s+") 
      e01_crc  = res[3]  ; e02_crc  = res[4]  ; e03_crc  = res[5]  ; e04_crc  = res[6]
      e05_crc  = res[7]  ; e06_crc  = res[8]  ; e07_crc  = res[9]  ; e08_crc  = res[10]
      e09_crc  = res[11] ; e10_crc  = res[12] ; e11_crc  = res[13] ; e12_crc  = res[14]



      aa = open(''+name[i]+'_crm_sci1.ap.mag.cal', 'w') 
      bb = open(''+name[i]+'_crm_sci1.ap.err.cal', 'w')
      cc = open(''+name[i]+'_crc_sci1.ap.mag.cal', 'w') 
      dd = open(''+name[i]+'_crc_sci1.ap.err.cal', 'w') 

      for j in range(len(xpsf)):
         idn    = j+1
         dis1   = np.sqrt( (xpsf[j]-x_crm)**2 + (ypsf[j]-y_crm)**2 ) 
         dis2   = np.sqrt( (xpsf[j]-x_crc)**2 + (ypsf[j]-y_crc)**2 )
         mat1   = (dis1 < 1.0) 
         mat2   = (dis2 < 1.0) 
         count1 = np.sum(mat1) 
         count2 = np.sum(mat2) 

         print(count1, name[i])
         if count1 == 1 and count2 == 1:
            mag1 = [m01_crm[mat1], m02_crm[mat1], m03_crm[mat1], m04_crm[mat1], 
                    m05_crm[mat1], m06_crm[mat1], m07_crm[mat1], m08_crm[mat1], 
                    m09_crm[mat1], m10_crm[mat1], m11_crm[mat1], m12_crm[mat1]]

            err1 = [e01_crm[mat1], e02_crm[mat1], e03_crm[mat1], e04_crm[mat1], 
                    e05_crm[mat1], e06_crm[mat1], e07_crm[mat1], e08_crm[mat1],
                    e09_crm[mat1], e10_crm[mat1], e11_crm[mat1], e12_crm[mat1]]

            mag2 = [m01_crc[mat2], m02_crc[mat2], m03_crc[mat2], m04_crc[mat2], 
                    m05_crc[mat2], m06_crc[mat2], m07_crc[mat2], m08_crc[mat2], 
                    m09_crc[mat2], m10_crc[mat2], m11_crc[mat2], m12_crc[mat2]]

            err2 = [e01_crc[mat2], e02_crc[mat2], e03_crc[mat2], e04_crc[mat2], 
                    e05_crc[mat2], e06_crc[mat2], e07_crc[mat2], e08_crc[mat2],
                    e09_crc[mat2], e10_crc[mat2], e11_crc[mat2], e12_crc[mat2]]

            mag1 = mag1 + (zero[i]-25.0) + 2.5*np.log10(expo[i]) 
            mag2 = mag2 + (zero[i]-25.0) + 2.5*np.log10(expo[i]) 


            x_crm_mat = np.float(x_crm[mat1])
            y_crm_mat = np.float(y_crm[mat1])
            if x_crm_mat > 20 and x_crm_mat < nelex-20 and y_crm_mat > 20 and y_crm_mat < neley-20 and np.float(mag1[4]) < 50:

               txt1 = fmt %(idn, xdrc[j], ydrc[j], x_crm[mat1], y_crm[mat1], 
                      mag1[0], mag1[1], mag1[2], mag1[3], mag1[4],  mag1[5], 
                      mag1[6], mag1[7], mag1[8], mag1[9], mag1[10], mag1[11])

               txt2 = fmt %(idn, xdrc[j], ydrc[j], x_crm[mat1], y_crm[mat1], 
                      err1[0], err1[1], err1[2], err1[3], err1[4],  err1[5],
                      err1[6], err1[7], err1[8], err1[9], err1[10], err1[11])

               txt3 = fmt %(idn, xdrc[j], ydrc[j], x_crc[mat2], y_crc[mat2], 
                      mag2[0], mag2[1], mag2[2], mag2[3], mag2[4],  mag2[5],
                      mag2[6], mag2[7], mag2[8], mag2[9], mag2[10], mag2[11])

               txt4 = fmt %(idn, xdrc[j], ydrc[j], x_crm[mat1], y_crm[mat1], 
                      err2[0], err2[1], err2[2], err2[3], err2[4],  err2[5],
                      err2[6], err2[7], err2[8], err2[9], err2[10], err2[11])

               aa.write(txt1)
               bb.write(txt2)
               cc.write(txt3)
               dd.write(txt4)
      aa.close()
      bb.close()
      cc.close()
      dd.close()






      if chan[i] == 'WFC' or chan[i] == 'UVIS':

         hdulist = fits.open(dir+''+name[i]+'_'+sfx[i]+'.fits')
         w = wcs.WCS(hdulist['sci',2].header, hdulist)
         hdulist.close()

         xpsf, ypsf = w.all_world2pix(ra, dec, 1)  

         res   = pd.read_csv(''+name[i]+'_crm_sci2.ap.mag', header=None, delimiter=r"\s+") 
         x_crm    = res[1]  ; y_crm    = res[2]
         m01_crm  = res[3]  ; m02_crm  = res[4]  ; m03_crm  = res[5]  ; m04_crm  = res[6]
         m05_crm  = res[7]  ; m06_crm  = res[8]  ; m07_crm  = res[9]  ; m08_crm  = res[10]
         m09_crm  = res[11] ; m10_crm  = res[12] ; m11_crm  = res[13] ; m12_crm  = res[14]

         res   = pd.read_csv(''+name[i]+'_crm_sci2.ap.err', header=None, delimiter=r"\s+") 
         e01_crm  = res[3]  ; e02_crm  = res[4]  ; e03_crm  = res[5]  ; e04_crm  = res[6]
         e05_crm  = res[7]  ; e06_crm  = res[8]  ; e07_crm  = res[9]  ; e08_crm  = res[10]
         e09_crm  = res[11] ; e10_crm  = res[12] ; e11_crm  = res[13] ; e12_crm  = res[14]

         res   = pd.read_csv(''+name[i]+'_crc_sci2.ap.mag', header=None, delimiter=r"\s+") 
         x_crc    = res[1]  ; y_crc    = res[2]
         m01_crc  = res[3]  ; m02_crc  = res[4]  ; m03_crc  = res[5]  ; m04_crc  = res[6]
         m05_crc  = res[7]  ; m06_crc  = res[8]  ; m07_crc  = res[9]  ; m08_crc  = res[10]
         m09_crc  = res[11] ; m10_crc  = res[12] ; m11_crc  = res[13] ; m12_crc  = res[14]

         res   = pd.read_csv(''+name[i]+'_crc_sci2.ap.err', header=None, delimiter=r"\s+") 
         e01_crc  = res[3]  ; e02_crc  = res[4]  ; e03_crc  = res[5]  ; e04_crc  = res[6]
         e05_crc  = res[7]  ; e06_crc  = res[8]  ; e07_crc  = res[9]  ; e08_crc  = res[10]
         e09_crc  = res[11] ; e10_crc  = res[12] ; e11_crc  = res[13] ; e12_crc  = res[14]

         aa = open(''+name[i]+'_crm_sci2.ap.mag.cal', 'w')
         bb = open(''+name[i]+'_crm_sci2.ap.err.cal', 'w')
         cc = open(''+name[i]+'_crc_sci2.ap.mag.cal', 'w')
         dd = open(''+name[i]+'_crc_sci2.ap.err.cal', 'w')

         for j in range(len(xpsf)):
            idn    = j+1
            dis1   = np.sqrt( (xpsf[j]-x_crm)**2 + (ypsf[j]-y_crm)**2 ) 
            dis2   = np.sqrt( (xpsf[j]-x_crc)**2 + (ypsf[j]-y_crc)**2 ) 
            mat1   = (dis1 < 1.0) 
            mat2   = (dis2 < 1.0) 
            count1 = np.sum(mat1) 
            count2 = np.sum(mat2) 

            if count1 == 1 and count2 == 1:
               mag1 = [m01_crm[mat1], m02_crm[mat1], m03_crm[mat1], m04_crm[mat1],
                       m05_crm[mat1], m06_crm[mat1], m07_crm[mat1], m08_crm[mat1],
                       m09_crm[mat1], m10_crm[mat1], m11_crm[mat1], m12_crm[mat1]]

               err1 = [e01_crm[mat1], e02_crm[mat1], e03_crm[mat1], e04_crm[mat1], 
                       e05_crm[mat1], e06_crm[mat1], e07_crm[mat1], e08_crm[mat1],
                       e09_crm[mat1], e10_crm[mat1], e11_crm[mat1], e12_crm[mat1]]

               mag2 = [m01_crc[mat2], m02_crc[mat2], m03_crc[mat2], m04_crc[mat2], 
                       m05_crc[mat2], m06_crc[mat2], m07_crc[mat2], m08_crc[mat2],
                       m09_crc[mat2], m10_crc[mat2], m11_crc[mat2], m12_crc[mat2]]

               err2 = [e01_crc[mat2], e02_crc[mat2], e03_crc[mat2], e04_crc[mat2], 
                       e05_crc[mat2], e06_crc[mat2], e07_crc[mat2], e08_crc[mat2],
                       e09_crc[mat2], e10_crc[mat2], e11_crc[mat2], e12_crc[mat2]]

               mag1 = mag1 + (zero[i]-25.0) + 2.5*np.log10(expo[i]) 
               mag2 = mag2 + (zero[i]-25.0) + 2.5*np.log10(expo[i]) 

               x_crm_mat = np.float(x_crm[mat1])
               y_crm_mat = np.float(y_crm[mat1])

               if x_crm_mat > 20 and x_crm_mat < nelex-20 and y_crm_mat > 20 and y_crm_mat < neley-20 and np.float(mag1[4]) < 50:

                  txt1 = fmt %(idn, xdrc[j], ydrc[j], x_crm[mat1], y_crm[mat1], 
                         mag1[0], mag1[1], mag1[2], mag1[3], mag1[4],  mag1[5],
                         mag1[6], mag1[7], mag1[8], mag1[9], mag1[10], mag1[11])

                  txt2 = fmt %(idn, xdrc[j], ydrc[j], x_crm[mat1], y_crm[mat1], 
                         err1[0], err1[1], err1[2], err1[3], err1[4],  err1[5],
                         err1[6], err1[7], err1[8], err1[9], err1[10], err1[11])

                  txt3 = fmt %(idn, xdrc[j], ydrc[j], x_crc[mat2], y_crc[mat2], 
                         mag2[0], mag2[1], mag2[2], mag2[3], mag2[4],  mag2[5],
                         mag2[6], mag2[7], mag2[8], mag2[9], mag2[10], mag2[11])

                  txt4 = fmt %(idn, xdrc[j], ydrc[j], x_crm[mat1], y_crm[mat1], 
                         err2[0], err2[1], err2[2], err2[3], err2[4],  err2[5],
                         err2[6], err2[7], err2[8], err2[9], err2[10], err2[11])

                  aa.write(txt1)
                  bb.write(txt2)
                  cc.write(txt3)
                  dd.write(txt4)
         aa.close()
         bb.close()
         cc.close()
         dd.close()











