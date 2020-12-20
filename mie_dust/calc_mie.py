import numpy as np
import matplotlib.pyplot as plt
import miepython
from matplotlib import gridspec
from astropy.io import fits, ascii
import glob
import pdb
import os
from joblib import Memory
from astropy.io import fits, ascii
from astropy.table import Table 
from scipy.interpolate import interp1d
import pkg_resources



labelDict = {"Mg2SiO4": "Mg$_2$SiO$_4$ (Forsterite)",
             "SiO2": "SiO$_2$ (quartz)",
             "Fe2SiO4":"Fe$_2$SiO$_4$ (fayalite)",
             "Fe":"Fe (Metallic Iron)",
             "MgSiO3":"MgSiO$_3$ (enstatite)",
             "Al2O3":"Al$_2$O$_3$ (corundum)",
             "C":"C (Carbon)"}
# Analytic Approximation
# def q_absorb(x, n_i, n_r):
#     num = 24*x*n_r*n_i
#     den = (n_r**2 + 2)**2
#     return num/den
#
# def q_scat(x,n_r):
#     num = 8*x**4*(n_r**2 -1)**2
#     den = 3*(n_r**2 +2)**2
#     return num/den
#
#
# def q_ext(x,n_i, n_r):
#     q_out = q_absorb(x,n_i,n_r) + q_scat(x,n_r)
#     return q_out
    
def q_ext_full(x,n_i, n_r):
    x_in = np.array(x)
    n_complex = np.array(n_r) - np.array(n_i) * 1j
    qext, qsca, qback, g = miepython.mie(n_complex,x_in)
    return qext

memory = Memory(cachedir='cache')
q_ext_mem = memory.cache(q_ext_full)

def all_opt_coeff_full(x,n_i,n_r):
    x_in = np.array(x)
    n_complex = np.array(n_r) - np.array(n_i) * 1j
    qext, qsca, qback, g = miepython.mie(n_complex,x_in)
    return qext, qsca, qback, g

all_opt_coeff_mem = memory.cache(all_opt_coeff_full)


def reshape_dotprod(mieResult,npoint,nwav,weights):
    mieResult2D = np.reshape(mieResult,(npoint,nwav))
    return np.dot(weights,mieResult2D)

def q_ext_lognorm_full(wavel,rad=1.0,n=complex(1.825,-1e-4),logNorm=False,
                  npoint=128,pdfThreshold=0.001,s=0.5):
    """
    Calculates the Mie extinction cross section Q_ext as a function of wavelength
    
    Arguments
    ------------------
    rad: float
        The radius of the particle
    wavel: float,arr
        The array of wavelengths to evaluate
    """
    sz = 2. * np.pi * rad/np.array(wavel)
    
    if logNorm == True:
        ## Generate an array of points along the distribution
        ## Size multiplier
        nwav = sz.size
        ## Size to evaluate lognormal weights
        ## Make a linear space from threshold to threshold in the PDF
        lowEval, highEval = invLognorm(s,rad,pdfThreshold)
        
        sizeEvalExtra = np.logspace(np.log(lowEval),np.log(highEval),base=np.e,num=(npoint+1))
        sizeEval = sizeEvalExtra[0:-1] ## extra point is for differentials
        dSize = np.diff(sizeEvalExtra)
        
        weights = lognorm(sizeEval,s,rad) * dSize
        sumWeights = np.sum(weights)
        if (sumWeights < 0.8) | (sumWeights > 1.1):
            print(r'!!!!!!Warning, PDF weights not properly sampling PDF!!!!')
            pdb.set_trace()
        weights = weights / sumWeights
        ## Arrange the array into 2D for multiplication of the grids
        sizeMult = (np.tile(sizeEval/rad,(nwav,1))).transpose()
        sizeArr = np.tile(sz,(npoint,1))
        
        eval2D = sizeMult * sizeArr
        finalEval = eval2D.ravel() ## Evaluate a 1D array
        
        n_2d = np.tile(n,(npoint,1))
        final_n = n_2d.ravel()
    else:
        finalEval = np.array(sz)
        final_n = n
    
    ## Make all points with sz > 100. the same as 100
    ## Otherwise the miescatter code quits without any warning or diagnostics
    highp = finalEval > 100.
    finalEval[highp] = 100.
    
    qext, qsca, qback, g = miepython.mie(final_n,finalEval)
    
    if logNorm == True:
        ## Now sum by the weights
        finalQext = reshape_dotprod(qext,npoint,nwav,weights)
        final_qscat = reshape_dotprod(qsca,npoint,nwav,weights)
        final_qback = reshape_dotprod(qback,npoint,nwav,weights)
        final_g = reshape_dotprod(g,npoint,nwav,weights)
    else:
        finalQext = qext
        
    return finalQext, final_qscat, final_qback, final_g

q_ext_lognorm = memory.cache(q_ext_lognorm_full)

def invLognorm(s,med,pdfThreshold):
    """ 
    Calculates the X values for a Log-normal distribution evaluated at specific PDF values
    Arguments
    ------------------
    s: float
        The sigma (scale value) of the log-normal distribution
    med: float
        The median particle size
    pdfThreshold: float
        The PDF threshold at which to find the x values
    """
    mu = np.log(med)
    z = np.log(s * np.sqrt(2. * np.pi) * pdfThreshold)
    sqrtPart = np.sqrt((2 * s**2 - 2 * mu)**2 - 4. * mu**2 - 8. * s**2 * z)
    lowEval = np.exp(mu - s**2 - 0.5 * sqrtPart)
    highEval = np.exp(mu - s**2 + 0.5 * sqrtPart)
    
    return lowEval, highEval

def lognorm(x,s,med):
    """
    Calculates a log-normal size distribution
    
    Arguments
    ------------------
    x: arr
        The input particle size
    s: float
        The sigma value
    med: float
        The median particle size
    """
    mu = np.log(med)
    y = 1. / (s* x * np.sqrt(2.*np.pi)) * np.exp(-0.5*((np.log(x)-mu)/s)**2)
    return y


def get_index_refrac(wav,material='Fe2SiO4'):
    """
    Return the index of refraction for a given wavelength and material
    
    Parameters
    ----------
    wav: float or numpy array
        Wavelength in microns to evaluate
    material: str
        Name of the material to look up
    """
    opt_data_path = 'optical_dat/{}[s].dat'.format(material)
    full_opt_data_path = pkg_resources.resource_filename('mie_dust',opt_data_path)
    
    dat = ascii.read(full_opt_data_path)
    
    f_k = interp1d(dat['wl(um)'],dat['k'])
    f_n = interp1d(dat['wl(um)'],dat['n'])
    
    k = f_k(wav)
    n = f_n(wav)
    
    return k,n

def get_mie_coeff(wav,r=0.1,material='Fe2SiO4'):
    """
    Return the Mie coefficients for a given radius and wavelength (single particle size)
    Assumes homogeneous spherical particles
    
    Parameters
    ----------
    wav: float or numpy array
        Wavelength in microns to evaluate
    r: float or numpy array
        Radii of the particles in microns
    material: str
        Name of the material to look up
    """
    
    x = 2. * np.pi * r/wav
    
    k, n = get_index_refrac(wav,material=material)
    
    qext, qsca, qback, g = all_opt_coeff_full(x, k, n)
    
    return qext, qsca, qback, g

def plot_all(r=0.1,plotType='extinction',distribution='single'):
    compList = glob.glob('optical_data/*.dat')

    fig, ax = plt.subplots(figsize=(5,4))

    np.random.seed(0)

    for ind,oneComp in enumerate(compList):
        dat = ascii.read(oneComp)
        thisLabel = (os.path.basename(oneComp)).split('[s]')[0]
        
        x = 2. * np.pi * r/dat['wl(um)'] 
        if distribution == 'single':
            q_ext_1 = q_ext(x, dat['k'], dat['n'])
        elif distribution == 'lognormal':
            complexN = np.array(dat['n']) - np.array(dat['k']) * 1j
            result = q_ext_lognorm(dat['wl(um)'],rad=r,n=complexN,logNorm=True,
                                    npoint=24)
            q_ext_1 = result[0]
        else:
            raise Exception("Unrecognized distribution")
        
        if plotType == 'extinction':
            yPlot = q_ext_1
        elif plotType == 'depth':
            pts = (dat['wl(um)'] > 0.45) & (dat['wl(um)'] < 0.9)  # Kepler
            
            if np.sum(pts) == 0:
                f_q_ext = interp1d(dat['wl(um)'],q_ext_1)
                normFactor = 0.5/ f_q_ext(0.675)
            else:
                normFactor = 0.5 / np.mean(q_ext_1[pts])
            depth = normFactor * q_ext_1
            yPlot = depth
            
            binWave, binErr = get_snr()
            f = interp1d(dat['wl(um)'],depth)
            bin_depth = f(binWave)
        else:
            raise Exception("Unrecognized plot type {}".format(plotType))
        
        if plotType == 'extinction':
            offset = ind * 0.8 + 2.
        else:
            offset = ind * 1.0
        
        ptsInfo, = ax.plot(dat['wl(um)'],yPlot + offset,label=labelDict[thisLabel])
        
        ax.text(3.0,offset + 0.5,labelDict[thisLabel],color=ptsInfo.get_color())
        
        if plotType == 'depth':
            simY = bin_depth + np.random.randn(len(bin_depth)) * binErr
            ax.errorbar(binWave,simY + offset,yerr=binErr,color=ptsInfo.get_color(),
                        fmt='o',capsize=10)
    
    ax.set_xlim(0., 13.)

    ax.set_xlabel("Wavelength ($\mu$m)")
    

    if distribution == 'single':
        thisTitle = "Single-radius R={} $\mu$m ".format(r)
    else:
        thisTitle = "Log-Normal R={} $\mu$m ".format(r)
    
    ax.set_title(thisTitle)
    
    #ax.legend(loc='upper left')
    rName = "{:.2f}".format(r).replace('.','p')
    if plotType == 'extinction':
        ax.set_ylim(0,12)
        ax.set_ylabel("Q$_{ext}$")
        outName = "extinction_{}_{}".format(rName,distribution)
    else:
        ax.set_ylim(-0.2,7.8)
        ax.set_ylabel("Depth (%)")
        outName = "depth_{}_{}".format(rName,distribution)
    
    fig.savefig('plots/png/{}.png'.format(outName),bbox_inches='tight')
    fig.savefig('plots/pdf/{}.pdf'.format(outName),bbox_inches='tight')
    plt.close(fig)
