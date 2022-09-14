function checkParameters(p)

    @unpack Nx,Ns,Nz,Xo,So,Pbo,Sbo = p

    # Check provided initial conditions 
    Nx == length(Xo)  || error("Number of Xo initial conditions should be ", Nx)
    Ns == length(So)  || error("Number of So initial conditions should be ", Ns)
    Nx == length(Pbo) || error("Number of Pbo initial conditions should be ", Nx)
    Ns == length(Sbo) || error("Number of Sbo initial conditions should be ", Ns)

    return
end