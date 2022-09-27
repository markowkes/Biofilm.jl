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
function analyzeBiofilm(sol,p,t; makePlot=false)

    @unpack Nx,Ns,Nz,XNames,SNames,Title = p

    println("Analyzing ",Title)

    # Print titles to REPL
    printBiofilmTitles(p)
    
    for tn in t
        # Unpack solution 
        X,S,Pb,Sb,Lf=unpack_solution(sol,p,tn)

        # Print values to REPL
        printBiofilmValues(tn,X,S,Pb,Sb,Lf,p)

        # Make plot of biofilm variables at this time
        if makePlot
            makeBiofilmPlots(tn,Pb,Sb,Lf,p)
        end
    end

    return nothing
end

function printBiofilmTitles(p)
    @unpack XNames,SNames = p
    # Build output string
    str=@sprintf(" %8s |"," Time  ")
    map( (x) -> str*=@sprintf(" %8.8s |",x),XNames)
    map( (x) -> str*=@sprintf(" %8.8s |",x),SNames)
    map( (x) -> str*=@sprintf(" min,max(%8.8s) |",x),XNames)
    map( (x) -> str*=@sprintf(" min,max(%8.8s) |",x),SNames)
    str*=@sprintf(" %8s"," Lf [μm] ")
    # Print string
    println(str)
    return
end

function printBiofilmValues(t,X,S,Pb,Sb,Lf,p)
    @unpack Nx,Ns = p
    # Build output string
    str=@sprintf(" %8.3f |",t)
    map( (x)   -> str*=@sprintf(" %8.3g |",x),X)
    map( (x)   -> str*=@sprintf(" %8.3g |",x),S)
    for i in 1:Nx
        map( (x,y) -> str*=@sprintf(" %8.3g,%8.3g |",x,y),minimum(Pb[i,:]),maximum(Pb[i,:]))
    end
    for i in 1:Ns
        map( (x,y) -> str*=@sprintf(" %8.3g,%8.3g |",x,y),minimum(Sb[i,:]),maximum(Sb[i,:]))
    end
    map( (x)   -> str*=@sprintf(" %8.3g",x),1e6*Lf)
    # Print string
    println(str)
    return
end