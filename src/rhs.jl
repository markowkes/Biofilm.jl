
function biofilmRHS!(dsol,sol,p_r,t) 
    # Process param and range inputs
    p=p_r[1]
    r=p_r[2]

    # Unpack structs 
    @unpack Nx,Ns,Nz,rho,Kdet,mu = p

    # Split sol into dependent variables
    X,S,Pb,Sb,Lf=sol[r.X],sol[r.S],sol[r.Pb],sol[r.Sb],sol[r.Lf]
    Lf=Lf[1] # Convert length 1 vector into float
    Pb=reshape(Pb,Nx,Nz)
    Sb=reshape(Sb,Ns,Nz)

    # Compute particulate concentration from volume fractions
    Xb=similar(Pb)
    for j=1:Nx
        Xb[j,:] = rho[j]*Pb[j,:]
    end

    # Update grid
    z=range(0.0,Lf,Nz+1)
    zm=0.5*(z[1:Nz]+z[2:Nz+1])
    dz=z[2]-z[1]
    g=biofilmGrid(z,zm,dz)

    # Compute intermediate variables 
    fluxS = computeFluxS(S,Sb,p,g)              # Flux of substrate in biofilm
    μb    = computeMu_biofilm(Sb,Xb,Lf,t,p,g)   # Growthrates in biofilm
    μt    = computeMu_tank(S,X,Lf,t,p,g)        # Growthrates in tank
    V     = computeVel(μb,Pb,p,g)               # Velocity of particulates
    Vdet  = Kdet*Lf^2                           # Detachment velocity
    fluxP = computeFluxP(Pb,V,Vdet,p)           # Flux of particulates in biofilm

    # Compute RHS's
    dsol[r.X] =dXdt(t,X,S,Xb,Lf,Vdet,μt,p,g)    # Tank particulates
    dsol[r.S] =dSdt(t,X,S,Lf,Pb,fluxS,μt,p)     # Tank substrates
    dsol[r.Pb]=dPbdt(μb,Sb,Pb,fluxP,p,g)        # Biofilm particulates 
    dsol[r.Sb]=dSbdt(μb,Xb,fluxS,p,g)           # Biofilm substrates
    dsol[r.Lf]=dLfdt(V,Vdet)                    # Biofilm thickness
    return nothing
end

# RHS of tank particulates
function dXdt(t,X,S,Xb,Lf,Vdet,μt,p,g)
    @unpack Nx,mu,Q,V,A,src = p
    dX=similar(X)
    for j in 1:Nx
        dX[j] = ( μt[j]*X[j]             # Growth
                - Q*X[j]/V               # Flow out
                + Vdet*A*Xb[j,end]/V     # From biofilm
                + src[j](S,X,p)[1] )     # Source term
    end
    return dX
end

# RHS of tank substrates
function dSdt(t,X,S,Lf,Pb,fluxS,μt,p) 
    @unpack Ns,Q,Sin,V,A,Yxs = p
    dS = zeros(Ns); 
    for k in 1:Ns                                 
        dS[k] = ( Q.*Sin[k](t)/V              # Flow in
                - Q.*S[k]     /V              # Flow out
                - A.*fluxS[k,end]/V           # Flux into biofilm
                - sum(μt.*X./Yxs[:,k]) ) # Used by growth
        #if p.neutralization == true
        #    dS(k) = dS(k) - p.kB(k)*p.kdis(k)*p.rho(1)*Pb(1,end); # Neutralization
        #end
    end
    return dS
end

# RHS of biofilm particulates 
function dPbdt(μb,Sb,Pb,fluxPb,p,g) 
    @unpack Nx,Nz,src,rho = p
    @unpack dz = g
    netFlux= (fluxPb[:,2:end]-fluxPb[:,1:end-1])/dz # Flux in/out
    growth = μb.*Pb                                 # Growth
    Source = zeros(Nx,Nz)
    println("Pb   ",size(Pb),typeof(Pb))
    println("Sb   ",size(Sb),typeof(Sb))
    println("rho   ",size(rho),typeof(rho))
    for j in 1:Nx
        for i in 1:Nz
            Source[j,i] = src[j](Sb[:,i],Pb[i]*rho[j],p)/rho[j]
        end
    end
    dPb  = growth - netFlux + Source;
    # Return RHS as a column vector
    dPb=reshape(dPb,Nx*Nz,1)
    return dPb
end

# RHS of biofilm substrates 
function dSbdt(μb,Xb,fluxS,p,g)
    @unpack Nx,Ns,Nz,src,rho,Yxs = p
    @unpack dz = g
    netFlux= (fluxS[:,2:end]-fluxS[:,1:end-1])/dz # Diffusion flux
    growth = zeros(Ns,Nz)
    for k in 1:Ns
        for j in 1:Nx
            growth[k,:] = growth[k,:] + μb[j,:].*Xb[j,:]./Yxs[j,k] # Used by growth
        end
    end
    dSb = netFlux - growth
    # Return RHS as a column vector
    dSb=reshape(dSb,Ns*Nz,1)
    return dSb
end

# RHS of biofilm thickness
function dLfdt(V,Vdet)
    Vfilm = V[end]          # Growth velocity at top of biofilm
    dLf = Vfilm - Vdet    # Surface Velocity 
    return [dLf]
end