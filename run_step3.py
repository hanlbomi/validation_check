def run_step3(dir, FID, med, wid, mbri, mfai, mhist, bs, cra, cdec, rad, cini, cfin, mini, mfin):

   import pandas as pd
   import numpy as np
   import matplotlib.pyplot as plt
   from astropy.io import fits
   import os
   import glob
   import astropy.wcs as wcs
   from astropy import units as u
   from astropy.coordinates import SkyCoord

   res  = pd.read_csv('../check_img.dat', header = None, delimiter=r"\s+", skiprows=3) # 
   fil     = res[1]

   uniq    = np.flip(np.unique(fil))
   if res[4][0] == 'IR':
      uniq     = np.unique(fil)
   dum      = (fil == uniq[FID-1])

   fil     = np.array(fil[dum])
   expo    = np.array(res[2][dum])
   inst    = np.array(res[3][dum])
   chan    = np.array(res[4][dum])
   name    = np.array(res[8][dum])
   newname = np.array(res[9][dum])
   zero    = np.array(res[12][dum])
   sfx     = np.array(res[13][dum])






   fig = plt.figure(figsize=(15./2.54, 15./2.54))
   ax = fig.add_axes([0.12, 0.10, 0.86, 0.88])

   xrini = 0.7 
   xrfin = 1.9
   yrini = 26.5
   yrfin = 18

   if cini != 'None':
      xrini = cini
   if cfin != 'None':
      xrfin = cfin
   if mini != 'None':
      yrini = mini
   if mfin != 'None':
      yrfin = mfin

   if chan[0] == 'IR':
      xrini = 0.7 
      xrfin = 2.2  
      yrini = 25.5    
      yrfin = 17   

      if cini != 'None':
         xrini = cini
      if cfin != 'None':
         xrfin = cfin
      if mini != 'None':
         yrini = mini
      if mfin != 'None':
         yrfin = mfin

   csize = 15

   ax.set_xlim(xrini, xrfin, auto=False)
   ax.set_ylim(yrini, yrfin)
   ax.tick_params(direction='in', length=4.0, labelsize=csize)
   ax.tick_params(which='minor',top=True,right=True,direction='in')
   ax.tick_params(axis='x', which='both', top='on')
   ax.tick_params(axis='y', which='both', right='on')
   ax.minorticks_on()
   ax.set_xlabel('C (r = 0.4 - 1.25 FWHM)',size=csize)
   ax.set_ylabel('F'+str(fil[0])+'W',size=csize)





   cflc = [] ; mflc = []
   xflc = [] ; yflc = []
   ra   = [] ; dec  = [] ; sep = []

   for i in range(len(name)):
      print(name[i], fil[i])

      res  = pd.read_csv(''+name[i]+'_crm_sci1.ap.mag', header = None, delimiter=r"\s+")
      xdum = res[1] ; ydum = res[2]
      m1   = res[4] ; m2   = res[5] 
      m3   = res[6] ; m4   = res[7]
      mdum = m4 + (zero[i]-25) + 2.5*np.log10(expo[i])


      hdulist = fits.open(dir+name[i]+'_'+sfx[i]+'.fits')
      w1 = wcs.WCS(hdulist[('sci',1)].header, hdulist)
      hdulist.close()
      tra, tdec = w1.all_pix2world(xdum, ydum, 1)
      c1   = SkyCoord(ra = tra*u.deg, dec = tdec*u.deg, frame='icrs')
      c2   = SkyCoord(ra = cra*u.deg, dec = cdec*u.deg, frame='icrs')
      tsep = c1.separation(c2)*u.deg    
      tsep = np.array(tsep)*60.0        

      cflc = np.append(cflc, m1-m2)
      mflc = np.append(mflc, mdum)
      xflc = np.append(xflc, xdum)
      yflc = np.append(yflc, ydum)
      ra   = np.append(ra,   tra)
      dec  = np.append(dec,  tdec)
      sep  = np.append(sep,  tsep)


      if chan[i] == 'WFC' or chan[i] == 'UVIS':
         res  = pd.read_csv(''+name[i]+'_crm_sci2.ap.mag', header = None, delimiter=r"\s+")
         xdum = res[1] ; ydum = res[2]
         m1   = res[4] ; m2   = res[5] 
         m3   = res[6] ; m4   = res[7] 
         mdum = m4 + (zero[i]-25) + 2.5*np.log10(expo[i])

         hdulist = fits.open(dir+name[i]+'_'+sfx[i]+'.fits')
         w1 = wcs.WCS(hdulist[('sci',2)].header, hdulist)
         hdulist.close()
         tra, tdec = w1.all_pix2world(xdum, ydum, 1)
         c1   = SkyCoord(ra = tra*u.deg, dec = tdec*u.deg, frame='icrs')
         c2   = SkyCoord(ra = cra*u.deg, dec = cdec*u.deg, frame='icrs')
         tsep = c1.separation(c2)*u.deg    
         tsep = np.array(tsep)*60.0         

         cflc = np.append(cflc, m1-m2)
         mflc = np.append(mflc, mdum)
         xflc = np.append(xflc, xdum)
         yflc = np.append(yflc, ydum)
         ra   = np.append(ra,   tra)
         dec  = np.append(dec,  tdec)
         sep  = np.append(sep,  tsep)


   var = (sep > rad) & (mflc < 40) 
   ax.scatter(cflc[var], mflc[var], marker = 'o', c = 'gray', s = 5, edgecolor = "none", alpha=0.5, zorder=1)

   from hist import histx
   from hist import histy
   xini = xrini-0.1
   yini = xrfin+0.1
   var  = (sep > rad) & (mflc < mhist) & (mflc > mbri) 
   xhi  = histx(xini, yini, bs)
   yhi  = histy(cflc[var], xini, yini, bs)

   dum  = (xhi > 1) & (xhi < 1.5)
   peak = np.max(yhi[dum])
   scail= (mhist -1 - yrfin)/peak
   ax.step(xhi, mhist-scail*yhi, where='mid', lw=1.5, c='r', zorder=2, alpha=0.5)

   ax.plot([med-wid, med-wid], [yrini, yrfin], c='blue', ls='--', alpha=0.7, zorder=3)
   ax.plot([med+wid, med+wid], [yrini, yrfin], c='blue', ls='--', alpha=0.7, zorder=3)
   xfill = [med-wid, med+wid, med+wid, med-wid]
   yfill = [mfai, mfai, mbri, mbri]
   ax.fill(xfill, yfill, 'deepskyblue', alpha=0.3, edgecolor="none", zorder=0)
   ax.text(0.6, 0.9, 'From FLC images', fontsize=csize, transform=ax.transAxes)

   
   psf = (sep > rad) & (mflc < mfai) & (mflc > mbri) & (cflc > med-wid) & (cflc < med+wid)
   ra  = np.array(ra[psf])
   dec = np.array(dec[psf])
   mflc= np.array(mflc[psf])
   hdulist = fits.open('../'+str(fil[0])+'_sci.fits') 
   w1 = wcs.WCS(hdulist[0].header, hdulist)
   hdulist.close()
   xd, yd = w1.all_world2pix(ra, dec, 1)

   xdrc = []
   ydrc = []
   map  = []
   for i in range(len(xd)):
      dis = np.sqrt( (xd[i]-xd)**2 + (yd[i]-yd)**2 )
      mat = (dis < 1.0)
      count = np.sum(mat)
      if xd[i] > 0:
         xdrc = np.append(xdrc,  xd[i])
         ydrc = np.append(ydrc,  yd[i])
         map  = np.append(map ,  mflc[i]) 
         xd[mat] = -10






   res1  = pd.read_csv('../Dolphot.dat',  header=None, delimiter=r"\s+", skiprows=3)
   xfit  = res1[5] ; yfit  = res1[6]
   if FID == 1:
      mfit = res1[8]
   if FID == 2:
      mfit = res1[14]

   xx_mat   = [] 
   yy_mat   = []
   mfit_mat = []
   for i in range(len(xdrc)):
      dis = np.sqrt( (xdrc[i]-xfit)**2 + (ydrc[i]-yfit)**2)
      mat = (dis < 1.0)
      count = np.sum(mat)
      if count == 1:
         xx_mat   = np.append(xx_mat,   xfit[mat])
         yy_mat   = np.append(yy_mat,   yfit[mat])
         mfit_mat = np.append(mfit_mat, mfit[mat])




   xx_passed   = []
   yy_passed   = []
   mfit_passed = []
   for i in range(len(xx_mat)):

      dis   = np.sqrt( (xx_mat[i]-xfit)**2 + (yy_mat[i]-yfit)**2)
      mat1  = (dis <= 30.) ; count1 = np.sum(mat1)
      mat2  = (dis <= 15.) ; count2 = np.sum(mat2)
      mat3  = (dis <= 10.) ; count3 = np.sum(mat3)
      mdum1 = np.array(mfit[mat1])
      mdum2 = np.array(mfit[mat2])
      mdum3 = np.array(mfit[mat3])
      ord1  = np.argsort(mdum1)
      ord2  = np.argsort(mdum2)
      ord3  = np.argsort(mdum3)
      mdum1 = mdum1[ord1]
      mdum2 = mdum2[ord2]
      mdum3 = mdum3[ord3]

      mdum11 = 999
      mdum22 = 999
      mdum33 = 999
      if count1 > 1.5:
         mdum11 = mdum1[1]
      if count2 > 1.5:
         mdum22 = mdum2[1]
      if count3 > 1.5:
         mdum33 = mdum3[1]
      if mfit_mat[i] == mdum1[0]:
         xx_passed   = np.append(xx_passed,   xx_mat[i])
         yy_passed   = np.append(yy_passed,   yy_mat[i])
         mfit_passed = np.append(mfit_passed, mfit_mat[i])



   from mk_reg import mk_reg
   mk_reg(''+str(fil[0])+'_temp.reg', xx_passed, yy_passed, 15, 1) 
   print('N of PSF candidates : ', len(xx_passed))
   print('')
   print('Check '+str(fil[0])+'c.png')
   print('')

   fig.savefig(''+str(fil[0])+'c.png', dpi=300)

