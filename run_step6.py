def run_step6(FID):

   import pandas as pd
   import numpy as np
   import os
   from matplotlib.backends.backend_pdf import PdfPages
   import matplotlib.pyplot as plt

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




   aa = open(''+str(fil[0])+'_cor_all.dat', 'w')
   aa.write('#ID  Name            X(drc)    Y(drc)    X(flc)     Y(flc)   m(crm)   e(crm)\n')
   aa.write('#                                                            r=5pix   r=5pix\n')


   with PdfPages(''+str(fil[0])+'_cor_all.pdf') as pdf:

      for i in range(len(fil)):

         for s in range(1, 3):
            res1  = pd.read_csv(''+str(name[i])+'_crm_sci'+str(s)+'.ap.mag.cal', header = None, delimiter=r"\s+")
            res2  = pd.read_csv(''+str(name[i])+'_crm_sci'+str(s)+'.ap.err.cal', header = None, delimiter=r"\s+")
            res3  = pd.read_csv(''+str(name[i])+'_crc_sci'+str(s)+'.ap.mag.cal', header = None, delimiter=r"\s+")
            res4  = pd.read_csv(''+str(name[i])+'_crc_sci'+str(s)+'.ap.err.cal', header = None, delimiter=r"\s+")
            idn    = res1[0]
            xdrc   = res1[1]  ; ydrc   = res1[2]
            xflc   = res1[3]  ; yflc   = res1[4]
            mcrm1  = res1[5]  ; ecrm1  = res2[5]  ;  mcrc1  = res3[5]  ; ecrc1  = res4[5] 
            mcrm2  = res1[6]  ; ecrm2  = res2[6]  ;  mcrc2  = res3[6]  ; ecrc2  = res4[6] 
            mcrm3  = res1[7]  ; ecrm3  = res2[7]  ;  mcrc3  = res3[7]  ; ecrc3  = res4[7] 
            mcrm4  = res1[8]  ; ecrm4  = res2[8]  ;  mcrc4  = res3[8]  ; ecrc4  = res4[8] 
            mcrm5  = res1[9]  ; ecrm5  = res2[9]  ;  mcrc5  = res3[9]  ; ecrc5  = res4[9] 
            mcrm6  = res1[10] ; ecrm6  = res2[10] ;  mcrc6  = res3[10] ; ecrc6  = res4[10] 
            mcrm7  = res1[11] ; ecrm7  = res2[11] ;  mcrc7  = res3[11] ; ecrc7  = res4[11] 
            mcrm8  = res1[12] ; ecrm8  = res2[12] ;  mcrc8  = res3[12] ; ecrc8  = res4[12] 
            mcrm9  = res1[13] ; ecrm9  = res2[13] ;  mcrc9  = res3[13] ; ecrc9  = res4[13] 
            mcrm10 = res1[14] ; ecrm10 = res2[14] ;  mcrc10 = res3[14] ; ecrc10 = res4[14] 
            mcrm11 = res1[15] ; ecrm11 = res2[15] ;  mcrc11 = res3[15] ; ecrc11 = res4[15] 
            mcrm12 = res1[16] ; ecrm12 = res2[16] ;  mcrc12 = res3[16] ; ecrc12 = res4[16] 


            xini = 0.13 
            xfin = 0.99
            yini = 0.11
            yfin = 0.99
            csize= 16
            crit = 1.0  

            rad = [0.05, 0.10, 0.15] 
            rad = np.array(rad)

            nele = len(xflc)
            for j in range(nele):
               print(name[i], idn[j])
               fig = plt.figure(figsize=(20./2.54, 15./2.54))
               ax  = fig.add_axes([xini, yini, xfin-xini, yfin-yini])

               mcrm = [mcrm1[j],mcrm2[j],mcrm3[j],mcrm4[j],mcrm5[j],mcrm6[j],mcrm7[j],mcrm8[j],mcrm9[j],mcrm10[j],mcrm11[j],mcrm12[j]]
               ecrm = [ecrm1[j],ecrm2[j],ecrm3[j],ecrm4[j],ecrm5[j],ecrm6[j],ecrm7[j],ecrm8[j],ecrm9[j],ecrm10[j],ecrm11[j],ecrm12[j]]
               mcrc = [mcrc1[j],mcrc2[j],mcrc3[j],mcrc4[j],mcrc5[j],mcrc6[j],mcrc7[j],mcrc8[j],mcrc9[j],mcrc10[j],mcrc11[j],mcrc12[j]]
               ecrc = [ecrc1[j],ecrc2[j],ecrc3[j],ecrc4[j],ecrc5[j],ecrc6[j],ecrc7[j],ecrc8[j],ecrc9[j],ecrc10[j],ecrc11[j],ecrc12[j]]

               mcrm = np.array(mcrm) #
               ecrm = np.array(ecrm) # 
               mcrc = np.array(mcrc) # 
               ecrc = np.array(ecrc) # 

               xrini = 0.0
               xrfin = 0.8
               yrini = np.max(mcrc)-0.5
               yrfin = np.min(mcrc)-0.3

               ax.set_xlim(xrini, xrfin)
               ax.set_ylim(yrini, yrfin)
               ax.xaxis.set_ticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
               ax.tick_params(direction='in', length=4.0, labelsize=csize) 
               ax.tick_params(axis='x', which='both', top='on')
               ax.tick_params(axis='y', which='both', right='on')
               ax.set_ylabel('Mag',size=csize)
               ax.set_xlabel('Aperture radius ["]',size=csize)
               ax.tick_params(which='minor',top=True,right=True,direction='in')
               ax.tick_params(which='minor',top=True,right=True,direction='in')
               ax.minorticks_on()

               ax.plot([xrini, xrfin], [mcrc[4], mcrc[4]], linestyle='--', c='silver', zorder=0)
               ax.plot([xrini, xrfin], [mcrc[9], mcrc[9]], linestyle='--', c='silver', zorder=0)
               ax.plot([rad[4], rad[4]], [yrini, yrfin], linestyle='--', c='silver', zorder=0)
               ax.plot([rad[9], rad[9]], [yrini, yrfin], linestyle='--', c='silver', zorder=0)

               ax.errorbar(rad, mcrc, yerr = ecrc, c='k', linestyle=' ', zorder=1)
               ax.plot(    rad, mcrc, c='b', zorder=2)
               ax.scatter( rad, mcrc, c='b', marker='o', s=15, zorder=2)

               dum = (mcrm < 40)
               ax.errorbar(rad[dum], mcrm[dum], yerr = ecrm[dum], c='k', linestyle=' ', zorder=3)
               ax.plot(    rad[dum], mcrm[dum], c='r', zorder=4)
               ax.scatter( rad[dum], mcrm[dum], c='r', marker='o', s=15, zorder=4)


               ax.text(0.08, 0.9 , 'ID=%3d' %(idn[j]), fontsize=csize, transform=ax.transAxes)
               ax.text(0.88, 0.9 , '%3d' %fil[i],  fontsize=csize, transform=ax.transAxes)
               name_all = name[i]+'_sci'+str(s)
               ax.text(0.70, 0.85, '%9s' %name_all, fontsize=csize, transform=ax.transAxes)
               ax.text(0.33, 0.45, 'XY$_{flc}$ = %9.3f, %9.3f' %(xflc[j], yflc[j]), 
               fontsize=csize, transform=ax.transAxes)
               ax.text(0.33, 0.38, 'XY$_{drc}$ = %9.3f, %9.3f' %(xdrc[j], ydrc[j]), 
               fontsize=csize, transform=ax.transAxes)
               ax.text(0.33, 0.31, 'm$_{crm}$ = %7.3f (5pix)' %mcrm[4], 
               fontsize=csize, transform=ax.transAxes, c='r')
               ax.text(0.33, 0.24, 'm$_{crc}$ = %7.3f (5pix), %7.3f (10pix)' %(mcrc[4], mcrc[9]), 
               fontsize=csize, transform=ax.transAxes, c='b')
               aa.write('%3d %15s %9.3f %9.3f %9.3f %9.3f %7.3f %7.3f     O\n' 
                         %(idn[j], name_all, xdrc[j], ydrc[j], xflc[j], yflc[j], mcrm[4], ecrm[4]) )

               pdf.savefig()
               plt.close()

            aa.write(' \n')

   aa.close()

   print( )
   os.system('scp '+str(fil[0])+'_cor_all.dat '+str(fil[0])+'_cor_all_vis.dat')
   print('Check '+str(fil[0])+'_cor.pdf')
   print('Edit '+str(fil[0])+'_cor_all_vis.dat')
   print( )

