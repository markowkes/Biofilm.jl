

using UnPack
using OrdinaryDiffEq
using DiffEqCallbacks
using Interpolations
 
 #Biofilm detachment function
function biofilmdetachment(integrator)

    sol=integrator.sol 
    p=integrator.p
    @unpack rho, Nx, Nz,outPeriod,plotPeriod= p
    @unpack rXt,rSt,rPb,rSb,rLf, A, V, prem = p

    # Convert solution to dependent variables
    t,Xt,St,Pb,Sb,Lf=unpack_solutionForPlot(sol,p)
    
    
    #set prem = percent biofilm remaining, find Lf initial and Lf remaining
    Lf = Lf[end]      
    Lfrem = Lf*prem

    state = "\n    Lf, Lfrem, amt detached"
    println(state)
    println(Lf)
    println(Lfrem)
    println(Lf-Lfrem)


    #compute grid and set variables
    grid = computeGrid(Lf, p)   #computeGrid only returns grid with initial values of Lf
    dz = grid.dz
    z = grid.z

    #find biofilm particulate mass before detachment
    state = "\n    before detachment, Xtank, Xbiofilm, sum of both"
    println(state)

    Xtank = Xt[:, end].*V       
    println(Xtank)
    Xbiofilm = sum(Pb, dims=2).*rho.*A.*dz 
    println(Xbiofilm)
    Xpredet = Xtank.+Xbiofilm
    println(Xpredet)
    
    #prepare biofilm remaining fractions, calculate fractions
    biof_rem = zeros(1, Nz)

    for n= 1:Nz
        frac_rem = max(0, min(1, (Lfrem - z[n])/dz))
        biof_rem[1, n] = frac_rem
    end

    #calculate amount of particulates dissolving, sum across rows, reshape for updating Xt

    Pbadd = (([1].-biof_rem).*Pb.*rho.*A.*dz./V)  
    Pbadd = sum(Pbadd, dims = 2)          
    Pbadd = reshape(Pbadd, :, 1)

    #update particulates in solution
    state = "\n    Pbadd, integrator.u[rXt], sum of both"
    println(state)
    
    println(Pbadd)
    println(integrator.u[rXt])
    integrator.u[rXt] = (integrator.u[rXt] + Pbadd)
    println(integrator.u[rXt])

    #update biofilm thickness
    integrator.u[rLf]= prem*integrator.u[rLf]  

    #calculate old and new locations of Nz points in grid
    xold = grid.zm
    xnew = (grid.zm).*prem
    newdz = dz.*prem

    #update particulate densities in remaining biofilm on new grid
        #println(Pb)
    interp_biofilm!(Pb,xold,xnew)
        #println(Pb)
    integrator.u[rPb] = Pb

    #update solution densities in remaining biofilm on new grid
    interp_biofilm!(Sb,xold,xnew)
    integrator.u[rSb] = Sb

    
    #find biofilm particulate mass after detachment
    state = "\n    after detachment, Xtank, Xbiofilm, sum of both"
    println(state)

    
    newPb = integrator.u[rPb]
    newPb=reshape(newPb,Nx,Nz)

    Xtank = integrator.u[rXt].*V     
    println(Xtank)
    Xbiofilm = sum(newPb, dims=2).*rho.*A.*newdz
    println(Xbiofilm)
    Xpostdet = Xtank.+Xbiofilm
    println(Xpostdet)


end
  
function interp_biofilm!(A,xold,xnew)

    # Loop over rows in A (each Pb or Sb)
    for n = 1:size(A,1)

        # Create interpolation operator
        interp_linear_extrap = linear_interpolation(xold, A[n,:], extrapolation_bc=Line())

        # Interpolate A[n,:] to xnew positions 
        A[n,:] = interp_linear_extrap(xnew)
    end
    return nothing
end

function detachmentTimes(tDeto, detPeriod, tFinal)
    #default value for tDeto is -1, results in no detachment events
    if tDeto == -1
        return []

    #default value for detPeriod is -1, if tDeto is provided but on detPeriod, one detachment event occurs at tDeto
    elseif detPeriod == -1
        return [tDeto]
    
    #if tDeto, detPeriod provided, then return range of values from tDeto to tFinal with detPeriod step
    else
        detTimes = range(start = tDeto, step = detPeriod, stop = tFinal)
        return detTimes
    end
end