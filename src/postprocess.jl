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
    p = checkInputs(p)
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
    p = checkInputs(p)
    @unpack Nx,Ns,Nz=p

    # Check t is a single time 
    length(t) > 1 && error("unpack_solution can only work with one time!")

    # Get solution at requested time 
    s=sol(t)

    # Compute ranges of dependent variables in sol array
    p = computeRanges(p)

    # Unpack solution
    Xt=s[p.rXt]
    St=s[p.rSt]
    Pb=s[p.rPb]
    Sb=s[p.rSb]
    Lf=s[p.rLf]

    # Reshape biofilm variables
    Pb=reshape(Pb,Nx,Nz)
    Sb=reshape(Sb,Ns,Nz)
    
    return Xt,St,Pb,Sb,Lf
end

"""
    biofilm_analyze(sol,p,times)
    biofilm_analyze(sol,p,times,makePlot=true)

Take solution from biofilm solver and outputs variabes and a plot of biofilm variables.

"""
function biofilm_analyze(sol,p,times; makePlot=false, plotSize=(1600,500))
    p = checkInputs(p)
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
            plt = biofilm_plot_film(sol([0,times[n]]),p)
            display(plt)
        end
    end

    return Xtout, Stout, Lfout
end

"""
    biofilm_movie(sol,p,times)
    biofilm_movie(sol,p,times,filename="anim.gif", fps=20)

Make a movie of the biofilm particulate volume fraction, solute concentration, and particulate growthrates at the specified times.

Optional arguments 
- filename: name and type of output, i.e., "biofilm.mp4", "biofilm.gif"
- framerate

Examples:

Create movie with t=0,1,...,10
```julia-repl
julia> biofilm_movie(sol,p,0:1:10)
```
Create movie with specified filename and framerate
```julia-repl
julia> biofilm_movie(sol,p,0:1:10,filename="biofilm.gif",fps=10)
```
"""
function biofilm_movie(sol,p,times; filename="anim.gif", fps=20)
    p = checkInputs(p)
    # Check times
    minimum(times) >= minimum(sol.t) || 
        error("Minimum time must be >= to ",minimum(sol.t))
    maximum(times) <= maximum(sol.t) || 
        error("Maximum time must be >= to ",maximum(sol.t))

    # Make animation
    anim = @animate for t in times 
        biofilm_analyze(sol,p,t,makePlot=true)
    end

    # Save annimation 
    gif(anim,filename,fps=fps)
end

"""
    biofilm_sol2csv(sol,p)
    biofilm_sol2csv(sol,p; filename)
Saves the solution of biofilm problem (sol) to a CSV file biofilm.csv 
The filename can be changed with the optional parameter

# Examples
```julia-repl
julia> biofilm_sol2csv(sol,p)
julia> biofilm_sol2csv(sol,p,filename="myfile.csv")
```
"""
function biofilm_sol2csv(sol,p; filename="biofilm.csv")
    p = checkInputs(p)
    @unpack Nx,Ns,Nz,XNames,SNames = p

    # Open file 
    open(filename, "w") do io
        
        # Write header 
        ##############
        write(io,"t,")
        # Tank Particulates
        for j = 1:Nx
            write(io, XNames[j]*"_t, ")
        end
        # Tank Solutes
        for j = 1:Ns
            write(io, SNames[j]*"_t, ")
        end
        # Biofilm Particulates
        for j = 1:Nx
            for i=1:Nz
                write(io, XNames[j]*"_b_"*string(i)*", ")
            end
        end
        # Biofilm Solutes
        for j = 1:Ns
            for i=1:Nz
                write(io, SNames[j]*"_b_"*string(i)*", ")
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
