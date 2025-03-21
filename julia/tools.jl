using Revise
using Healpix
using LinearAlgebra
using StaticArrays
using Base.Threads
using WignerD
using BenchmarkTools
using NPZ
using Plots
using PyPlot
using PyCall
using DataStructures
using HDF5
using Falcons

hp = pyimport("healpy")


function truncate_alm(alm; lmax_in, lmax_out, mmax_out)

    """
    This function trancates the alm.

    Args:
        alm: healpy's assumed alm file

        lmax_in: alm's lmax
    
        lmax_out: Output's lmax. It should be lower than lmax_in.
    
        mmax_out: Output's mmax. It should be lower than lmax_out.

    Returns:
        alm_new: Truncated alm.
    """
    @show size_out = alm_idx(lmax_out,mmax_out,lmax_out)
    alm_new = zeros(ComplexF64, 3, size_out)
    for m in 0:mmax_out
        idx_start_out = alm_idx(m,m,lmax_out)
        idx_stop_out = alm_idx(lmax_out, m,lmax_out)
        idx_start_in = alm_idx(m, m,lmax_in)
        idx_stop_in = alm_idx(lmax_out, m, lmax_in)
        
        alm_new[1, idx_start_out:idx_stop_out] = alm[1, idx_start_in:idx_stop_in]
        alm_new[2, idx_start_out:idx_stop_out] = alm[2, idx_start_in:idx_stop_in]
        alm_new[3, idx_start_out:idx_stop_out] = alm[3, idx_start_in:idx_stop_in]
    end
    
    return alm_new
end

function alm_idx(l, m::Integer, lmax::Integer)
    return Int(m * (2 * lmax + 1 - m) // 2 + l)+1
end