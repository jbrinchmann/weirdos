import numpy as np
from astropy.table import Table
from astropy.io import fits
import os.path
from scipy.sparse import csr_matrix

# Check whether resizing is available
try:
    from skimage.transform import resize
    resizing_available = True
except:
    resizing_available = False


# The location of the data
ROOTDIR = '/data2/jarle/MUSE/Weirdoes/'

def _get_suffix(bin, VoronoiSN):
    """Simple helper routine to create filenames suffixes"""
    
    if VoronoiSN is not None:
        suffix = "SN{0:03d}".format(VoronoiSN)
    else:
        suffix = "bin{0:02d}".format(bin)

    return suffix

def _voronoi_filename(name, VoronoiSN):
    """Helper routine to get the name of the file with the results of Voronoi binning"""
    suffix = _get_suffix(bin, VoronoiSN)
    voronoi_file = os.path.join(ROOTDIR, name,'voronoi-results-'+suffix+'.fits')

    return voronoi_file
    


def load_info():
    """Load the information on the weirdos

    Returns
    -------
    info : dict
        The dict has two keys: Name and Filename
           Name: The name of the galaxy
       Filename: The name of the Datacube belonging to this galaxy.
    """

    infofile = os.path.join(ROOTDIR, 'weirdos.info')
    info = Table().read(infofile, format='ascii')

    return info



def read_fits_columns(fname, what):
    """Read one or more columns from a fits binary table

    Parameters
    ----------
    fname : str
            The filename of the binary FITS table.
     what : str or list of str
            The quantities to load from the table.

    Returns
    -------
    res : dict
          Each requested column in the FITS table is returned
          in this dict with the column name as key

    Examples
    --------
    >>> wanted = ['H_ALPHA_FLUX']
    >>> d = read_fits_columns(<filename>, wanted)
    >>> print(d['H_ALPHA_FLUX'])
    """
    
    res = {}
    what_list = what if isinstance(what, list) else [what]
    try:
        hdul = fits.open(fname, ignore_missing_end=True)
        t = hdul[1].data # For platefit the first extension contains the data

        
        for w in what_list:
            res[w] = t.field(w)

        hdul.close()
    except:
        print("Something went wrong in reading the columns")
        raise


    return res

def load_autoz_results(name, bin=3, VoronoiSN=None,
                           toget=['ZMAP', 'ZMAP_FULL', 'VERR',
                                      'VERR_FULL', 'PROB', 'PROB_FULL']):
    """Load the results of running AutoZ

    AutoZ was run for the galaxies using both spatial binning and Voronoi binning.
    In the former case a binning of 3 was used, while in the latter a target
    S/N of 250 was used. This routine will load the Voronoi binned results
    if VoronoiSN is not set to None, otherwise the binned results.


    Parameters
    ----------
    name : str
           The name of the galaxy whose results you want.
     bin : int
           The spatial binning requested. Default: 3
    VoronoiSN : int
           The target S/N used for Voronoi binning. Only 250 is currently done.
           Default: None
    toget : list of str
           The quantities to lead. The default is to get ZMAP (the redshift map),
           ZMAP_FULL (like ZMAP but at original size if binned), VERR (the
           velocity uncertainty), VERR_FULL (c.f. ZMAP_FULL), PROB (an assessment
           of the probability of the cross-correlation quality),
           PROB_FULL (c.f. ZMAP_FULL)


    Returns
    -------
    result : dict
             A dict with the results for each element in toget, with toget used
             as keys.


    Example
    -------
    >>> name = 'SDSS004533-010606'
    >>> a = load_autoz_results(name, bin=3)
    >>> a.keys()                                                                                   
      dict_keys(['ZMAP', 'ZMAP_FULL', 'VERR', 'VERR_FULL', 'PROB', 'PROB_FULL'])
    >>> a['ZMAP_FULL'].shape                                                                       
      (316, 315)

    """

    suffix = _get_suffix(bin, VoronoiSN)
    fname = os.path.join(ROOTDIR, name, name+'-autoz-'+suffix+'.fits')

    hdul = fits.open(fname)
        
    result = dict()
    for t in toget:
        try:
            result[t] = hdul[t].data
        except:
            print("I did not find {0} .. skipping".format(t))
            pass

    return result



def load_voronoi_binning(name, VoronoiSN):
    """Load some of the information from Voronoi binning

    Parameters
    ----------
    name : str
           The name of the galaxy whose results you want.
     bin : int
           The spatial binning requested. Default: 3
    VoronoiSN : int
           The target S/N used for Voronoi binning. Only 250 is currently done.
           Default: None

    Returns
    -------
    result : dict
            This dictionary contains the main results of the 
            Voronoi binning and has the following keys:

            binnum: This is an image with same size as the original 
                    image (datacube spatial extent) and where the 
                    pixel value is the bin that pixel belongs to.
           snlimit: The target S/N used for Voronoi binning.
            header: The header (with WCS) belonging to binnum.
              mask: The mask used for the binning.
             n_bin: The number of bins
           reverse: A list where each entry contains the indices of
                    the pixels of each Voronoi bin belong to the 
                    image.
              area: The area of each Voronoi bin in pixels


    Example
    -------
    >>> name = 'SDSS004533-010606'
    >>> v = w.load_voronoi_binning(name, 250)
    
    """
    
    voronoi_file = _voronoi_filename(name, VoronoiSN)
    try:
        hdul = fits.open(voronoi_file)
        mainheader = hdul[0].header
        snlimit = mainheader['SNLIMIT']
        binnum = hdul['BINNUMBER'].data
        hdr = hdul['BINNUMBER'].header
        mask = hdul['MASK'].data
        hdul.close()
        n_bin = np.max(binnum).astype(int)

        # And we will need the reverse lookup. This is inefficient
        # in python compared to IDL.
        r = []
        area = []
        for i in range(n_bin+1):
            inds = np.where(binnum == i)
            r.append(inds)
            area.append(len(inds))
        

        result = {'binnum': binnum, 'snlimit': snlimit, 'header': hdr,
                    'mask': mask, 'n_bin': n_bin, 'reverse': r, 'area': area}
    except:
        print("An error occurred when trying to load the Voronoi binning")
        raise
        result = None

    return result




def voronoi_to_image(v, value):
    """Given a 1D array, use the Voronoi binning structure to convert
    this to a 2D image

    The issue here is when you run Voronoi binning you get a single vector,
    not an image, where the value corresponds to what was had from the 
    Voronoi binning. To make an image you need to reapply the Voronoi binning.
    
    Parameters
    ----------
         v : dict returned from load_voronoi_binning
             This contains the Voronoi binning information.
     value : ndarray
             The quantity to reform into an image.

    Returns
    -------
    im : ndarray
         An image of the quantity value.

    Example
    -------

    This example uses multiple functions here

    >>> name = 'SDSS004533-010606'
    >>> v = w.load_voronoi_binning(name, 250)
    >>> fname = 'SDSS004533-010606/Platefit/platefit-results-Voronoi-SN250.fit'
    >>> d = read_fits_columns(fname, ['H_ALPHA_EQW'])
    >>> im = voronoi_to_image(v, d)

    """

    # First create the image.
    im = np.zeros_like(v['binnum'])
    
    # We skip the first bin
    for i in range(1, v['n_bin']+1):
        if (v['area'][i] > 0):
            im[v['reverse'][i]] = value[i-1]

    return im

        
def load_platefit_results(name, what, bin=3, VoronoiSN=None, full_size=False):
    """Load the results of a platefit run

    Parameters
    ----------
       name : str
              The name of the galaxy
       what : str or list of strings
              The quantity(ies) to load from the platefit results file
        bin : int
              The spatial binning requested. Default: 3
  VoronoiSN : int
              The target S/N used for Voronoi binning. Only 250 is currently done.
              Default: None
  full_size : bool
              If true, resize the image of the platefit result (in the binned case)
              to match the full size of the datacube.

    Returns
    -------
    res : dict
          Each requested platefit result is returned
          in this dict with the column name as key

    Example
    -------

    Create OIII/Hb map of one galaxy

    >>> name = 'SDSS004533-010606'         
    >>> need = ['H_BETA_FLUX', 'OIII_5007_FLUX']
    >>> d = w.load_platefit_results(name, need, bin=3, full_size=True) 
    >>> o3_hb = np.log10(d['OIII_5007_FLUX']/d['H_BETA_FLUX'])
    >>> plt.imshow(o3_hb, vmin=-0.5, vmax=0.6)
    """


    if (full_size) and (not resizing_available):
        print("To be able to resize you need to install scikit-image")
        full_size=False
    
    suffix = _get_suffix(bin, VoronoiSN)
    if (VoronoiSN is not None):
        suffix = 'Voronoi-'+suffix
    fname = os.path.join(ROOTDIR, name, 'Platefit', 'platefit-results-'+suffix+'.fit')
    
    res = read_fits_columns(fname, what)

    # A challenge here is that the shape of the fits files. So we load
    # the mask for this object and reshapes the results to this.
    maskfile = os.path.join(ROOTDIR, name, 'mask.fits')
    hdul = fits.open(maskfile)
    nx = hdul[0].header['NAXIS1']
    ny = hdul[0].header['NAXIS2']
    hdul.close()
    
    # But the data might be binned?
    if (VoronoiSN is None):
        nx_bin = int(np.ceil(nx/bin))
        ny_bin = int(np.ceil(ny/bin))
        for k in res.keys():
            x = np.reshape(res[k], (ny_bin, nx_bin))
            if full_size:
                # Interpolate to full size.
                x = resize(x, (ny, nx), anti_aliasing=False, mode='constant')
            res[k] = x

    else:
        # Now we have a Voronoi binned dataset. This is a bit more involved because
        # now we need to get the Voronoi binning information.
        
        v = load_voronoi_binning(name, VoronoiSN)
        if (v is None):
            print("There was a problem with loading the binning.")
            print("1D data will be returned which you need to interpret yourself")
        else:
            for k in res.keys():
                x = voronoi_to_image(v, res[k])
                res[k] = x
                
        
    return res


    
