# Convert solution (1D vec) into more meaningful dependent variables
# ## This method is optimized for plotting ##
# t,X,S,Lf = f(t), Pb,Sb = f(z)
function unpack_solutionForPlot(sol,p,r)
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

"""
    unpack_solution(sol,p,t)

Convert `sol` into named dependent variables at the givien time t.

# Example
```julia-repl 
julia> X_here,S_here,Pb_here,Sb_here,Lf_here=unpack_solution(sol,p,3)
```
"""
function unpack_solution(sol,p,t)
    @unpack Nx,Ns,Nz=p

    # Check t is a single time 
    length(t) > 1 && error("unpack_solution can only work with one time!")

    # Get solution at requested time 
    s=sol(t)

    # Compute ranges of dependent variables in sol array
    # sol=[X,S,Pb,S,Lf]
    nVar=0; 
    N=Nx;    rX =nVar+1:nVar+N; nVar+=N # X =u[rX]
    N=Ns;    rS =nVar+1:nVar+N; nVar+=N # S =u[rS]
    N=Nx*Nz; rPb=nVar+1:nVar+N; nVar+=N # Pb=u[rPb]
    N=Ns*Nz; rSb=nVar+1:nVar+N; nVar+=N # Sb=u[rSb]
    N=1;     rLf=nVar+1:nVar+N          # Lf=u[rLf]

    # Unpack solution
    X=s[rX]
    S=s[rS]
    Pb=s[rPb]
    Sb=s[rSb]
    Lf=s[rLf]

    # Reshape biofilm variables
    Pb=reshape(Pb,Nx,Nz)
    Sb=reshape(Sb,Ns,Nz)
    
    return X,S,Pb,Sb,Lf
end

"""
    analyzeBiofilm(sol,p,t)
    analyzeBiofilm(sol,p,t,makePlot=true)

Take solution from biofilm solver and outputs variabes and a plot of biofilm variables.

"""
function analyzeBiofilm(sol,p::param,t;makePlot::Bool=false)

    @unpack Nx,Ns,Nz,XNames,SNames,Title = p

    println("Analyzing ",Title)
    
    for tn in t
        # Unpack solution 
        X,S,Pb,Sb,Lf=unpack_solution(sol,p,tn)

        # Print to REPL
        printBiofilmSolution(t,X,S,Pb,Sb,Lf,p)

        # Make plot of biofilm variables 
        if makePlot
            # Biofilm grid
            z=range(0.0,Lf[end],Nz+1)
            zm=0.5*(z[1:Nz]+z[2:Nz+1])

            # Plot solution in biofilm at this time 
            makeBiofilmPlots(tn,zm,Pb,Sb,p)
        end
    end

    return nothing
end

"""
    printBiofilmSolution(t,X,S,Pb,Sb,Lf,p)

Print biofilm variables at time t to REPL

"""
function printBiofilmSolution(t,X,S,Pb,Sb,Lf,p)
    # Unpack params
    @unpack Nx,Ns,XNames,SNames,Title = p

    # Print solution values at this time
    printstyled("                     t = ",t,"                     \n",
    bold=true,underline=true)
    printstyled("Tank:\n",bold=true)
    for i in 1:Nx 
        println(@sprintf(" %20s           = %10f [g/m³]",XNames[i],X[i]))
    end
    for i in 1:Ns
        println(@sprintf(" %20s           = %10f [g/m³]",SNames[i],S[i]))
    end
    printstyled("Biofilm:\n",bold=true)
    for i in 1:Nx 
        println(@sprintf(" %20s (min,max) = %10f,%10f [g/m³]",XNames[i],minimum(Pb[i,:]),maximum(Pb[i,:])))
    end
    for i in 1:Ns
        println(@sprintf(" %20s (min,max) = %10f,%10f [g/m³]",SNames[i],minimum(Sb[i,:]),maximum(Sb[i,:])))
    end
    println(@sprintf(" %20s           = %10f [μm]","Lf",1e6*Lf[1]))
    return nothing 
end