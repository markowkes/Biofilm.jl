using Accessors 

"""
    computeGrid(Lf,p)

Creates a grid for biofilm of thickness Lf
"""
function computeGrid(Lf,p)
    @unpack Nz = p
    z=range(0.0,Lf,Nz+1)
    zm=0.5*(z[1:Nz]+z[2:Nz+1])
    dz=z[2]-z[1]
    return biofilmGrid(z,zm,dz)
end

"""
    t,X,S,Lf = f(t), Pb,Sb = f(z)
Convert solution (1D vec) into more meaningful dependent variables
## This method is optimized for plotting ##

"""
function unpack_solutionForPlot(sol,p)
    @unpack Nx,Ns,Nz=p
    @unpack rXt,rSt,rPb,rSb,rLf = p
    t=sol.t
    Xt=sol[rXt,:]
    St=sol[rSt,:]
    Pb=sol[rPb,end]
    Sb=sol[rSb,end]
    Lf=sol[rLf,:]

    # Reshape biofilm variables
    Pb=reshape(Pb,Nx,Nz)
    Sb=reshape(Sb,Ns,Nz)
    
    return t,Xt,St,Pb,Sb,Lf
end

""" 
    computeRanges(p)
Compute the location of dependent variables in solution vector 
"""
function computeRanges(p)
    @unpack Nx,Ns,Nz = p 
    nVar=0; 
    N=Nx;    rXt=nVar+1:nVar+N; nVar+=N # Xt=u[rXt]
    N=Ns;    rSt=nVar+1:nVar+N; nVar+=N # St=u[rSt]
    N=Nx*Nz; rPb=nVar+1:nVar+N; nVar+=N # Pb=u[rPb]
    N=Ns*Nz; rSb=nVar+1:nVar+N; nVar+=N # Sb=u[rSb]
    N=1;     rLf=nVar+1:nVar+N          # Lf=u[rLf]

    # Add ranges to parameter struct to be used in RHS calc
    @reset p[:rXt] = rXt
    @reset p[:rSt] = rSt
    @reset p[:rPb] = rPb
    @reset p[:rSb] = rSb
    @reset p[:rLf] = rLf
    return p
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

function printBiofilmValues(t,Xt,St,Pb,Sb,Lf,p)
    @unpack Nx,Ns = p
    # Build output string
    str=@sprintf(" %8.3f |",t)
    map( (x)   -> str*=@sprintf(" %8.3g |",x),Xt)
    map( (x)   -> str*=@sprintf(" %8.3g |",x),St)
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


