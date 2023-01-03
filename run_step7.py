def run_step7():

   import pandas as pd
   import numpy as np
   import os
   from matplotlib.backends.backend_pdf import PdfPages
   import matplotlib.pyplot as plt

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

   fil     = np.flip(np.unique(fil)) 


   res   = pd.read_csv('../olot.dat',  header=None, delimiter=r"\s+", skiprows=3)
   xfit  = res[5] ; yfit  = res[6]
   mi    = res[7] ; mi_re = res[8] ; mei   = res[9]
   mv    = res[13]; mv_re = res[14]; mev   = res[15]

   for i in range(len(fil)):

      delt = -2.5*np.log10(0.853) 

      aa = open(''+str(fil[i])+'_cor.dat', 'w')
      aa.write('#ID   Xdrc      Ydrc       N   Map(mean) Map(med) Map       stdev    error   Mfit      corr\n')
      aa.write('#                              (r=5pix) (r=5pix) (r=inf)                               (map - mfit)\n')

      res  = pd.read_csv(''+str(fil[i])+'_cor_all.dat', header = None, delimiter=r"\s+", skiprows=2)
      dum  = (res[8] == 'O')  
      idn  = np.array(res[0][dum])
      xdrc = np.array(res[2][dum])
      ydrc = np.array(res[3][dum])
      xflc = np.array(res[4][dum])
      yflc = np.array(res[5][dum])
      m5   = np.array(res[6][dum])
      me5  = np.array(res[7][dum])
      flag = np.array(res[8][dum])
      idn_uniq = np.unique(idn)

      print(idn)

      if fil[i] == fil1:
         mfit = mi_re    
      if fil[i] == fil2: 
         mfit = mv_re    

      for j in range(len(idn_uniq)):
         dum = (idn == idn_uniq[j])
         print(idn_uniq[j], np.sum(dum), np.array(m5[dum]))
         xdrc_uniq = np.mean(np.array(xdrc[dum]))
         ydrc_uniq = np.mean(np.array(ydrc[dum]))
         m5_mean   = np.mean(np.array(m5[dum]))
         m5_median = np.median(np.array(m5[dum]))
         me5_mean  = np.mean(np.array(me5[dum])) 
         m5_stdev  = np.std(np.array(m5[dum]))
         count     = np.sum(dum)
         minf      = m5_mean - delt

         dis = np.sqrt( (xdrc_uniq-xfit)**2 + (ydrc_uniq-yfit)**2 )
         mat = (dis < 1.0)
         corr = minf - mfit[mat] 
 
         if np.sum(mat) == 1:
            aa.write('%4d %9.3f %9.3f %3d %8.3f %9.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n' 
            %(idn_uniq[j], xdrc_uniq, ydrc_uniq, count, m5_mean, m5_median, minf, m5_stdev, me5_mean, mfit[mat], corr) )

      aa.close()
