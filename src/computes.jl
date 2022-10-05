# Fluxes of substrate due to diffusion: F=De*dSb/dz
function computeFluxS(S,Sb,p,g)
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
function computeMu_biofilm(Sb,Xb,Lf,t,p,g)
    @unpack Nx,Nz,mu = p 
    @unpack zm = g
    μb=zeros(Nx,Nz)
    for j in 1:Nx
        μb[j,:]=mu[j](Sb,Xb,Lf,t,zm,p)
    end
    return μb
end

# Growthrate for each particulate in tank
function computeMu_tank(S,X,Lf,t,p,g)
    @unpack Nx,Nz,mu = p 
    @unpack z = g
    μt=zeros(Nx)
    for j in 1:Nx
        μt[j]=mu[j](S,X,Lf,t,z[end],p)[1]
    end
    return μt
end

# Velocity due to growth in biofilm
function computeVel(μb,Sb,Pb,t,p,g)
    @unpack Nx,Nz,Ptot,src,rho = p
    @unpack dz = g
    # Velocities on faces of cells
    V=zeros(Nz+1) 
    # Start with zero velocity at wall -> integrate through the biofilm
    for i in 1:Nz
        V[i+1]=V[i]
        for j in 1:Nx
            # Add growth of particulates in this cell to velocity
            V[i+1] += μb[j,i].*Pb[j,i]*dz/Ptot
            # Add source of particulates in this cell to velocity
            V[i+1] += src[j](Sb[:,i],Pb[:,i]*rho[j],t,p)[1]/rho[j]*dz/Ptot
        end
    end
    
    return V
end

# Fluxes of particulate due to advection: F=V*phi;
function computeFluxP(Pb,V,Vdet,p)
    @unpack Nx,Nz = p
    # Fluxes
    # Bottom boundary - no flux condition -> nothing to do
    fluxP = zeros(Nx,Nz+1) # Fluxes on faces of cells
    for i in 2:Nz+1  # Interior and top faces
            fluxP[:,i]= V[i]*Pb[:,i-1] # V*phi_face (upwinded)
        end
    return fluxP
end
