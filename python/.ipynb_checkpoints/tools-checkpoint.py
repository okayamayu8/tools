import numpy as np
import healpy as hp

def truncate_alm(alm, lmax_in, lmax_out):
    """
    This function trancates the alm.

    Args:
        alm: healpy's assumed alm file

        lmax_in: alm's lmax
        
        lmax_out: Output's lmax. It should be lower than lmax_in.

    Returns:
        alm_new: Truncated alm.
    """
    size_out = hp.Alm.getsize(lmax_out)
    alm_new  = np.zeros((3,size_out), dtype=alm.dtype)
    for m in range(lmax_out+1):
        idx_start_out = hp.Alm.getidx(lmax_out, m, m)
        idx_stop_out  = hp.Alm.getidx(lmax_out, lmax_out, m)
        idx_start_in  = hp.Alm.getidx(lmax_in, m, m)
        idx_stop_in   = hp.Alm.getidx(lmax_in, lmax_out, m)
        alm_new[0,idx_start_out:idx_stop_out+1] = alm[0,idx_start_in:idx_stop_in+1]
        alm_new[1,idx_start_out:idx_stop_out+1] = alm[1,idx_start_in:idx_stop_in+1]
        alm_new[2,idx_start_out:idx_stop_out+1] = alm[2,idx_start_in:idx_stop_in+1]
    return alm_new

def gaussian_sidelobe(nside, fwhm, peak):
    sigma = fwhm/(2*np.sqrt(2*np.log(2)))
    npix = hp.nside2npix(nside)
    maps = np.zeros((3,npix), np.float64)
    pixs = np.linspace(0,npix-1,npix)
    theta,phi=hp.pix2ang(nside, np.int32(pixs))
    result = 1/(2*np.pi*sigma**2)*np.exp(-(theta - peak)**2/sigma**2)
    maps[0] = result
    maps[1] = result*np.cos(2*phi)
    maps[2] = -result*np.sin(2*phi)
    return maps