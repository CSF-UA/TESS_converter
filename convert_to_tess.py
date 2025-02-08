import os
import matplotlib
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import math
import time
import scipy.optimize as spo



def get_list_of_stars(sector):
    LOS = []
    f = open('tesscurl_sector_'+ str(sector) +'_lc.sh','r')
    S = f.readlines()[1:]
    f.close()
    for i in range(len(S)):
        a = S[i]
        a = S[i].split('-')
        LOS.append(int(a[6]))
    return LOS



def save_data(JD, mag, star, date, sector):
    
    X = [[JD[j], mag[j]] for j in range(len(JD))]
    X = np.array(X)
    start = path_C+'/'+str(sector)+'/'
    
    fig = plt.figure(1)
    fig.set_size_inches(15, 8)
    plt.plot(JD, mag, '.k', markersize = 2)
    plt.gca().invert_yaxis()
    plt.title('TIC '+str(star) + ' sector ' + str(sector), fontsize = 20)
    plt.xlabel('JD - 2 457 000', fontsize = 16)
    plt.ylabel('magnitude, mmag', fontsize = 16)
    matplotlib.rc('xtick', labelsize=16)
    matplotlib.rc('ytick', labelsize=16)
    plt.savefig(start+'TIC_'+str(star)+'_'+date+'_mag.png', dpi=100)
    plt.close()
    
    fname = start+'TIC_'+str(star)+'_'+date+'_mag.tess'
    np.savetxt(fname, X)



def clear_up(JD, mag):
    JD = list(JD)
    mag = list(mag)
    i = 0
    while i<len(JD):
        if math.isnan(mag[i]) or math.isnan(JD[i]):
            del mag[i]
            del JD[i]
        else:
            i += 1
    return np.array(JD), np.array(mag)



def convert_sector(sector):
    LOS = get_list_of_stars(sector)
    f = open('obs_start.txt', 'r')
    S = f.readlines()
    f.close()
    a = S[sector-1].split()
    date = a[0]
    x1 = a[1]
    x2 = a[2]
    
    for i in range(len(LOS)):
        try:
            sc = str(sector)
            sc = 's'+'0' * (4-len(sc)) + sc
            number = str(LOS[i])
            number = '0' * (16-len(number)) + number
            name = x1 +'-'+ sc +'-'+ number +'-0'+ x2 +'-s_lc.fits'
            fits_file = path_C+'/'+str(sector)+'/'+name
            fits.info(fits_file)
            fits.getdata(fits_file, ext=1).columns
            with fits.open(fits_file, mode="readonly") as hdulist:
                JD = hdulist[1].data['TIME']
                flux = hdulist[1].data['PDCSAP_FLUX']
            
            JD = list(JD)
            flux = list(flux)
            tf = [math.isnan(flux[j]) for j in range (len(flux))]
            j = 0
            while j<len(flux):
                if tf[j] == True:
                    del tf[j]
                    del flux[j]
                    del JD[j]
                else:
                    j += 1
            
            JD, flux = np.array(JD), np.array(flux)
            mag = -2.5*np.log10(flux)
            mag -= np.average(mag)
            mag *= 1000
            
            save_data(JD, mag, LOS[i], date, sector)
            
        except:
            print('ERROR')



path_C = os.getcwd()
convert_sector(79)