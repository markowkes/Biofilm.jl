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
function computeVel(μb,Pb,p,g)
    @unpack Nz,Ptot = p
    @unpack dz = g
    # Velocities on faces of cells
    V=zeros(Nz+1) 
    # Start with zero velocity at wall -> integrate through the biofilm
    for i in 1:Nz
        # Add growth of particulates in this cell to velocity
        V[i+1]=V[i] + sum(μb[:,i].*Pb[:,i]*dz/Ptot);
    end
    return V
end

# Fluxes of particulate due to advection: F=V*phi;
function computeFluxP(Pb,V,Vdet,p)
    @unpack Nx,Nz = p
    # Fluxes
    fluxP = zeros(Nx,Nz+1) # Fluxes on faces of cells
    for i in 2:Nz  # Interior faces
        fluxP[:,i]= V[i]*(Pb[:,i-1]+Pb[:,i])/2 # V*phi_face
    end
    # Bottom boundary - no flux condition -> nothing to do
    # Top boundary - use phi in top cell
    fluxP[:,Nz+1] = V[Nz+1]*(Pb[:,Nz]) #- Vdet*Pb(:,end);
    return fluxP
end

# Convert solution (1D vec) into more meaningful dependent variables
# t,X,S,Lf = f(t), Pb,Sb = f(z)
function unpack_solution(sol,p,r)
    @unpack Nx,Ns,Nz=p
    t=sol.t
    X=sol[r.X,:]
    S=sol[r.S,:]
    Pb=sol[r.Pb,end]
    Sb=sol[r.Sb,end]
    Lf=sol[r.Lf,:]

    # Reshape biofilm variables
    Pb=reshape(Pb,Nx,Nz)
    Sb=reshape(Sb,Ns,Nz)
    
    return t,X,S,Pb,Sb,Lf
end

# Return greatest common divisor of floats
function Base.gcd(a::Float64,b::Float64)

    # Check order of numbers
    a < b && return gcd(b, a)

    # Check for Inf 
    a ≈ Inf && return b
    
    # base case
    if abs(b) < eps(Float64)
        return a
    else
        return (gcd(b, a - floor(a / b) * b))
    
    end
end
