

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
    t,Xt,St,Pb,Sb,Lf=unpack_solutionForPlot(sol,p)      #rename this
    
    #set prem = percent biofilm remaining, find Lf initial and Lf remaining
    Lf = Lf[end]      
    Lfrem = Lf*prem

    #compute grid and set variables
    grid = computeGrid(Lf, p) 
    dz = grid.dz
    z = grid.z

    #prepare biofilm remaining fractions, calculate fractions
    biof_rem = zeros(1, Nz)

    for n= 1:Nz
        frac_rem = max(0, min(1, (Lfrem - z[n])/dz))
        biof_rem[1, n] = frac_rem
    end

    #update particulates in solution
    Pbadd = calcDissolve(Pb, dz, V, A, rho, biof_rem)
    integrator.u[rXt] = (integrator.u[rXt] + Pbadd)

    #calculate and update solutes added to solution
    #rho = 1 for Sb as rho is a particulate density parameter
    Sbadd = calcDissolve(Sb, dz, V, A, 1, biof_rem)      
    integrator.u[rSt] = (integrator.u[rSt] + Sbadd)

    #update biofilm thickness
    integrator.u[rLf]= prem*integrator.u[rLf]  

    #calculate old and new locations of Nz points in grid
    xold = grid.zm
    xnew = (grid.zm).*prem

    #update particulate and solution densities in remaining biofilm on new grid
    interp_biofilm!(Pb,xold,xnew)
    integrator.u[rPb] = Pb

    interp_biofilm!(Sb,xold,xnew)
    integrator.u[rSb] = Sb

    return nothing
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

function detachmentTimes(tDeto, detPeriod, tFinal, deviation)
    #default value for tDeto is -1, results in no detachment events
    if tDeto == -1
        return []

    #default value for detPeriod is -1, if tDeto is provided but not detPeriod, one detachment event occurs at tDeto
    elseif detPeriod == -1
        return [tDeto]
    
    #default value for deviation is -1, if tDeto and detPeriod provided but not deviation then return range of values from tDeto to tFinal with detPeriod step
    
    elseif deviation == -1
        detTimes = range(start = tDeto, step = detPeriod, stop = tFinal)
        return detTimes
    
    #if all are provided, return range of detachment times with +/- deviation on each value
    else
        #create times
        detTimes = range(start = tDeto, step = detPeriod, stop = tFinal)

        #change lazily defined vector returned by range function into an array
        detTimes = collect(detTimes)

        #create uniformly distributed random number -1 to 1 for each time in detTimes and multiply by deviation
        numberTimes = length(detTimes)
        randomchange = (rand(numberTimes).-0.5).*2*deviation

        #add random deviations to original times
        detTimes = detTimes .+ randomchange
        return detTimes
    
    end
end

function calcDissolve(Array, dz, V, A, rho, biof_rem)
    #calculate amount of particulates or solutes detaching, sum across rows, reshape for updating Xt or St
    
    ArrayAdd = (([1].-biof_rem).*Array.*rho.*A.*dz./V)  
    ArrayAdd = sum(ArrayAdd, dims = 2)          
    ArrayAdd = reshape(ArrayAdd, :, 1)

    return ArrayAdd
end 