# Fluxes of substrate due to diffusion: F=De*dSb/dz
@timeit to function computeFluxS(S,Sb,p,g)
    @unpack Ns,Nz,De,Daq,LL = p
    @unpack dz = g
    fluxS = zeros(Ns,Nz+1); # Fluxes on faces of cells
    for i in 2:Nz  # Interior faces
        fluxS[:,i]= De[:].*(Sb[:,i]-Sb[:,i-1])/dz;
    end
    # Bottom boundary - no flux condition -> nothing to do
    # Top boundary - flux matching between biofilm and boundary layer 
    S_top= ( (Daq*(dz/2).*S+De*LL.*Sb[:,Nz])
            ./ (Daq*(dz/2)+De*LL) ) 
    fluxS[:,Nz+1] .= De.*(S_top-Sb[:,Nz])/(dz/2)
    return fluxS
end

# Growthrate for each particulate in biofilm
@timeit to function computeMu_biofilm(Sb,Xb,Lf,t,p,g)
    @unpack Nx,Nz,mu = p 
    @unpack zm = g
    μb=zeros(Nx,Nz)
    for j in 1:Nx
        μb[j,:]=mu[j](Sb,Xb,Lf,t,zm,p)
    end
    return μb
end

# Growthrate for each particulate in tank
@timeit to function computeMu_tank(S,X,Lf,t,p,g)
    @unpack Nx,Nz,mu = p 
    @unpack z = g
    μt=zeros(Nx)
    for j in 1:Nx
        μt[j]=mu[j](S,X,Lf,t,z[end],p)[1]
    end
    return μt
end

# Velocity due to growth in biofilm
@timeit to function computeVel(μb,Pb,srcb,p,g)
    @unpack Nx,Nz,Ptot = p
    @unpack dz = g
    # Velocities on faces of cells
    @timeit to "zeros" V=zeros(Nz+1) 
    # Start with zero velocity at wall -> integrate through the biofilm
    for i in 1:Nz
        V[i+1]=V[i]
        #for j in 1:Nx
            # Add growth of particulates in this cell to velocity
            @timeit to "growth" V[i+1] += sum(μb[:,i].*Pb[:,i])*dz/Ptot
            # Add source of particulates in this cell to velocity
            @timeit to "source" V[i+1] += sum(srcb[:,i])*dz/Ptot
        #end
    end
    
    return V
end

# Fluxes of particulate due to advection: F=V*phi;
@timeit to function computeFluxP(Pb,V,Vdet,p)
    @unpack Nx,Nz = p
    # Fluxes
    # Bottom boundary - no flux condition -> nothing to do
    fluxP = zeros(Nx,Nz+1) # Fluxes on faces of cells
    for i in 2:Nz+1  # Interior and top faces
        fluxP[:,i]= V[i]*Pb[:,i-1] # V*phi_face (upwinded)
    end
    return fluxP
end

# Source within biofilm
@timeit to function computeSource_biofilm(Sb,Pb,t,p,g)
    @unpack Nx,Nz,Ptot,src,rho = p
    @unpack dz = g

    srcb = zeros(Nx,Nz)
    for j in 1:Nx
        for i in 1:Nz
            srcb[j,i] = src[j](Sb[:,i],Pb[:,i]*rho[j],t,p)[1]/rho[j]
        end
    end
    return srcb
end
