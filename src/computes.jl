"""
    computeS_top(St,Sb,p,g)

Compute solute concentration at top of biofilm
using flux matching between biofilm and boundary layer
"""
function computeS_top(St,Sb,p,g)
    @unpack Nz,Db,Dt,LL = p
    @unpack dz = g
    S_top= ( (Dt*(dz/2).*St+Db*LL.*Sb[:,Nz])
            ./ (Dt*(dz/2)+Db*LL) )
    return S_top 
end

# Fluxes of solute due to diffusion: F=De*dSb/dz
function computeFluxS(St,Sb,p,g)
    @unpack Ns,Nz,Db = p
    @unpack dz = g
    fluxS = zeros(Ns,Nz+1); # Fluxes on faces of cells
    for i in 2:Nz  # Interior faces
        fluxS[:,i]= Db[:].*(Sb[:,i]-Sb[:,i-1])/dz;
    end
    # Bottom boundary - no flux condition -> nothing to do
    # Top boundary - flux matching between biofilm and boundary layer 
    S_top = computeS_top(St,Sb,p,g)
    fluxS[:,Nz+1] .= Db.*(S_top-Sb[:,Nz])/(dz/2)
    return fluxS
end

# Growthrate for each particulate in biofilm
function computeMu_biofilm(Sb,Xb,Lf,t,p,g)
    @unpack Nx,Nz,mu = p 
    @unpack zm = g
    μb=zeros(Nx,Nz)
    for j in 1:Nx
        for i in 1:Nz
            μb[j,i]=mu[j](Sb[:,i],Xb[:,i],Lf,t,zm[i],p)
        end
    end
    return μb
end

# Growthrate for each particulate in tank
function computeMu_tank(St,Xt,Lf,t,p,g)
    @unpack Nx,Nz,mu = p 
    @unpack z = g
    μt=zeros(Nx)
    for j in 1:Nx
        μt[j]=mu[j](St,Xt,Lf,t,z[end],p)[1]
    end
    return μt
end

# Velocity due to growth in biofilm
function computeVel(μb,Sb,Pb,Lf,t,p,g)
    @unpack Nx,Nz,Ptot,srcX,rho = p
    @unpack zm,dz = g
    # Velocities on faces of cells
    V=zeros(Nz+1) 
    # Start with zero velocity at wall -> integrate through the biofilm
    for i in 1:Nz
        V[i+1]=V[i]
        for j in 1:Nx
            # Add growth of particulates in this cell to velocity
            V[i+1] += μb[j,i].*Pb[j,i]*dz/Ptot
            # Add source of particulates in this cell to velocity
            V[i+1] += srcX[j](Sb[:,i],Pb[:,i]*rho[j],Lf,t,zm[i],p)[1]/rho[j]*dz/Ptot
        end
    end
    
    return V
end

# Fluxes of particulate due to advection: F=V*phi;
function computeFluxP(Pb,V,p)
    @unpack Nx,Nz = p
    # Fluxes
    # Bottom boundary - no flux condition -> nothing to do
    fluxP = zeros(Nx,Nz+1) # Fluxes on faces of cells
    for i in 2:Nz+1  # Interior and top faces
        fluxP[:,i]= V[i]*Pb[:,i-1] # V*phi_face (upwinded)
    end
    return fluxP
end
