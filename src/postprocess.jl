""" 
    MeanBiofilmVarsWithTime(sol,p)

Get mean (throughout biofilm) of S̄b[Ns,Nt] and X̄b[Ns,Nt] as function of time
Inputs: sol at end of simulaton and parameters p

# Example
```julia-repl 
julia> Sb_t,Pb_t=MeanBiofilmVarsWithTime(sol,p)

```
"""
function MeanBiofilmVarsWithTime(sol,p)
    @unpack Nx,Ns,Nz=p
    times = sol.t; Nt=length(times)
    Sb_t=zeros(Ns,Nt)
    Pb_t=zeros(Nx,Nt)
    for i in eachindex(times)
        t=times[i]
        Xt,St,Pb,Sb,Lf = unpack_solution(sol,p,t)
        for k=1:Nz
            for j=1:Ns 
                Sb_t[j,i] += Sb[j,k]/Nz
            end
            for j=1:Nx 
                Pb_t[j,i] += Pb[j,k]/Nz
            end
        end
    end
    return Sb_t,Pb_t
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
    N=Nx;    rXt=nVar+1:nVar+N; nVar+=N # X =u[rXt]
    N=Ns;    rSt=nVar+1:nVar+N; nVar+=N # S =u[rSt]
    N=Nx*Nz; rPb=nVar+1:nVar+N; nVar+=N # Pb=u[rPb]
    N=Ns*Nz; rSb=nVar+1:nVar+N; nVar+=N # Sb=u[rSb]
    N=1;     rLf=nVar+1:nVar+N          # Lf=u[rLf]

    # Unpack solution
    Xt=s[rXt]
    St=s[rSt]
    Pb=s[rPb]
    Sb=s[rSb]
    Lf=s[rLf]

    # Reshape biofilm variables
    Pb=reshape(Pb,Nx,Nz)
    Sb=reshape(Sb,Ns,Nz)
    
    return Xt,St,Pb,Sb,Lf
end

"""
    analyzeBiofilm(sol,p,times)
    analyzeBiofilm(sol,p,times,makePlot=true)

Take solution from biofilm solver and outputs variabes and a plot of biofilm variables.

"""
function analyzeBiofilm(sol,p,times; makePlot=false, plotSize=(1600,500))

    @unpack Nx,Ns,Nz,XNames,SNames,Title = p

    println("Analyzing ",Title)

    # Preallocate arrays for output 
    Xtout =  Array{Float64}(undef,p.Nx,length(times))
    Stout =  Array{Float64}(undef,p.Ns,length(times))
    Lfout = Vector{Float64}(undef,     length(times))

    # Print titles to REPL
    printBiofilmTitles(p)
    
    for n in eachindex(times)
        # Unpack solution 
        Xt,St,Pb,Sb,Lf=unpack_solution(sol,p,times[n])

        # Store in output arrays
        Xtout[:,n] .= Xt[:]
        Stout[:,n] .= St[:]
        Lfout[  n]  = Lf[1] 

        # Print values to REPL
        printBiofilmValues(times[n],Xt,St,Pb,Sb,Lf,p)

        # Make plot of biofilm variables at this time
        if makePlot
            makeBiofilmPlots(times[n],Pb,Sb,Lf,p,plotSize)
        end
    end

    return Xtout, Stout, Lfout
end

"""
    movieBiofilm(sol,p,times)
    movieBiofilm(sol,p,times,filename="anim.gif", fps=20)

Make a movie of the biofilm particulate volume fraction, substrate concentration, and particulate growthrates at the specified times.

Optional arguments 
- filename: name and type of output, i.e., "biofilm.mp4", "biofilm.gif"
- framerate

Examples:

Create movie with t=0,1,...,10
```julia-repl
julia> movieBiofilm(sol,p,0:1:10)
```
Create movie with specified filename and framerate
```julia-repl
julia> movieBiofilm(sol,p,0:1:10,filename="biofilm.gif",fps=10)
```
"""
function movieBiofilm(sol,p,times; filename="anim.gif", fps=20)

    # Check times
    minimum(times) >= minimum(sol.t) || 
        error("Minimum time must be >= to ",minimum(sol.t))
    maximum(times) <= maximum(sol.t) || 
        error("Maximum time must be >= to ",maximum(sol.t))

    # Make animation
    anim = @animate for t in times 
        analyzeBiofilm(sol,p,t,makePlot=true)
    end

    # Save annimation 
    gif(anim,filename,fps=fps)
end

"""
    sol2csv(sol,filename,p)
Saves the solution of biofilm problem (sol) to a CSV file with name filename.
"""
function sol2csv(sol,filename,p)
    @unpack Nx,Ns,Nz,XNames,SNames = p

    # Open file 
    open(filename, "w") do io
        
        # Write header 
        ##############
        write(io,"t,")
        # Tank Particulates
        for j = 1:Nx
            write(io, XNames[j]*", ")
        end
        # Tank Substrates
        for j = 1:Ns
            write(io, SNames[j]*", ")
        end
        # Biofilm Particulates
        for j = 1:Nx
            for i=1:Nz
                write(io, XNames[j]*"_"*string(i)*", ")
            end
        end
        # Biofilm Substrates
        for j = 1:Ns
            for i=1:Nz
                write(io, SNames[j]*"_"*string(i)*", ")
            end
        end
        # Biofilm thickness 
        write(io,"Lf \n")

        # Write values 
        ##############
        for n in 1:length(sol.t)
            write(io, string(sol.t[n])*", ")
            for i in 1:(Nx + Ns + Nx*Nz + Ns*Nz + 1)
                write(io,string(sol.u[n][i])*", ")
            end
            write(io, "\n")
        end
    end;

end
