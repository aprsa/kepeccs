import numpy as np
import urllib
import json
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Get all face value widths and separation:
KIC, p0 = np.loadtxt('kepEBs.csv', usecols=(0, 1), delimiter=',', unpack=True)
KIC = KIC.astype(int)
P0 = dict(zip(KIC, p0))

def linseg(phi, a, b, c):
    retval = np.empty_like(phi)
    retval[phi <= phi0] = a*phi[phi <= phi0] + b
    retval[phi  > phi0] = c*(phi[phi > phi0]-phi0) + a*phi0 + b
    return retval

def sqseg(phi, a, b, c, d, e):
    retval = np.empty_like(phi)
    phi_left = phi[phi <= phi0]
    phi_right = phi[phi > phi0]
    retval[phi <= phi0] = a*phi_left**2 + b*phi_left + c
    retval[phi  > phi0] = d*phi_right**2 + e*phi_right + (a-d)*phi0**2 + (b-e)*phi0 + c
    return retval

def linsqseg(phi, a, b, c, d):
    retval = np.empty_like(phi)
    phi_left = phi[phi <= phi0]
    phi_right = phi[phi > phi0]
    retval[phi <= phi0] = a*phi_left + b
    retval[phi  > phi0] = c*(phi_right**2-phi0**2) + d*phi_right + (a-d)*phi0 + b
    if 2*c*phi0+d >= 0:
        retval += 1e4*(2*c*phi0+d)
    return retval

# pull data from the server:
resp = urllib.urlopen("http://keplerebs.villanova.edu/api/mark_stats")
data = json.loads(resp.read())
# kics = [d['kic'] for d in data]
# sanjay = np.loadtxt('sanjay.list', dtype=int)
# for s in sanjay:
#     if s not in kics:
#         print('KIC %08d is not in.' % s)
# exit()

for i, d in enumerate(data):
    kic = d['kic']
    try:
        P0[kic]
    except:
        continue
    
    # read in the data:
    try:
        ph, fl = np.loadtxt('kepLCs/%08d.00.lc.data' % kic, usecols=(1, 6), unpack=True)
    except:
        print('skipping KIC %d' % kic)
        continue

    ecl1start, ecl1center, ecl1end, ecl2start, ecl2center, ecl2end = d['ecl_marks']

    # we'll use each of these consecutively to try and fit a good model to them.

    # ingress of the first eclipse.
    if ecl1start > -1.5: # value -2 means no eclipse
        chi2 = []
        sds = []
        prod = []

        phi0 = ecl1start
        phi0s = np.linspace(phi0-(ecl1end-ecl1start)/5, phi0+(ecl1end-ecl1start)/5, 101)

        dphis = np.linspace(0.05, 0.45, 101)*(ecl1end-ecl1start)
        dphi = 0.2*(ecl1end-ecl1start)

        slopes = np.ones((len(phi0s), len(dphis)))*np.nan

        for pi, phi0 in enumerate(phi0s):
            for di, dphi in enumerate(dphis):
                # filter the data to the subinterval of interest:
                filt = (ph > phi0-dphi) & (ph < phi0+dphi)
                if len(ph[filt]) < 5:
                    continue

                popt, pcov = curve_fit(linsqseg, ph[filt], fl[filt])
                # chi2.append( ((fl[filt]-linsqseg(ph[filt], *popt))**2).sum() )
                # sds.append(2*popt[2]*phi0+popt[3]-popt[0])
                # prod.append(chi2[-1]*sds[-1])
                slopes[pi, di] = 2*popt[2]*phi0+popt[3]-popt[0]

            # minchi2idx = np.array(chi2).argmin()
            # minsdsidx = np.array(sds).argmin()
            # minprodidx = np.array(prod).argmin()

            # print('%08d %12.8f %12.8f %12.8f' % (kic, dphi, ecl1start, phi0s[minprodidx]))

        plt.imshow(slopes.T)
        plt.show()

        continue
        # exit()

        PLOT = False
        PLOT_TO_FILE = True
        if PLOT or PLOT_TO_FILE:
            # phi0 = phi0s[minchi2idx]
            # filt = (ph > phi0-dphi) & (ph < phi0+dphi)
            # popt, pcov = curve_fit(linsqseg, ph[filt], fl[filt])

            # plt.plot(ph[filt], fl[filt], 'bo')
            # plt.plot(np.sort(ph[filt]), linsqseg(np.sort(ph[filt]), *popt), 'g-', lw=2)
            # plt.plot([phi0,], linsqseg(np.array([phi0,]), *popt), 'gs', ms=10)

            # phi0 = phi0s[minsdsidx]
            # filt = (ph > phi0-dphi) & (ph < phi0+dphi)
            # popt, pcov = curve_fit(linsqseg, ph[filt], fl[filt])

            # plt.plot(ph[filt], fl[filt], 'bo')
            # plt.plot(np.sort(ph[filt]), linsqseg(np.sort(ph[filt]), *popt), 'r-', lw=2)
            # plt.plot([phi0,], linsqseg(np.array([phi0,]), *popt), 'rs', ms=10)

            phi0 = phi0s[minprodidx]
            filt = (ph > phi0-dphi) & (ph < phi0+dphi)
            popt, pcov = curve_fit(linsqseg, ph[filt], fl[filt])

            plt.plot(ph[filt], fl[filt], 'bo')
            plt.plot(np.sort(ph[filt]), linsqseg(np.sort(ph[filt]), *popt), 'r-', lw=2)
            plt.plot([phi0,], linsqseg(np.array([phi0,]), *popt), 'rs', ms=10)

            plt.axvline(ecl1start)
            plt.savefig('ie_figs/%08d.e1s.png' % kic)
            if PLOT:
                plt.show()
            else:
                plt.clf()

            # plt.plot(phi0s, chi2/np.array(chi2).max(), 'b-')
            # plt.plot(phi0s, sds/np.array(sds).max(), 'r-')
            # plt.show()
