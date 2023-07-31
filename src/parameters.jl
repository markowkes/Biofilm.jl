using Printf
using Accessors

"""
    p = checkInputs(p)

Runs checks on inputs and sets default values 
"""
function checkInputs(p)
        p = checkTypes_setDefs(p; verbose=false)
        checkParameters(p)
    return p
end

"""
    printInputs(d)

Print named tuple d 
"""
function printInputs(d)
    for key in collect(keys(d))
        @printf(" %-19s => %s \n",key,(d[key]))
    end
    return nothing
 end

 
"""
    checkTypes_setDefs(d)

Checks the parameters in the named tupple d for common errors
and sets default values if not already specified. 
"""
function checkTypes_setDefs(d; verbose=false)

    # Check values of each field 
    # !!! set to default if missing and default exists !!!
    err = false
    d,err = checkType_setDef(err,d,Float64,              :tFinal                                     )
    d,err = checkType_setDef(err,d,Float64,              :outPeriod                                  )
    d,err = checkType_setDef(err,d,Float64,              :tol                                        )
    d,err = checkType_setDef(err,d,Float64,              :plotPeriod,         default=d[:outPeriod]  )
    d,err = checkType_setDef(err,d,Bool,                 :makePlots,          default=true           )
    d,err = checkType_setDef(err,d,Float64,              :discontinuityPeriod,default=Inf            )
    d,err = checkType_setDef(err,d,String,               :optionalPlot,       default="growthrate"   )
    d,err = checkType_setDef(err,d,Tuple{Int64, Int64},  :plotSize,           default=(1600,1000)    )
    d,err = checkType_setDef(err,d,String,               :Title                                      )
    d,err = checkType_setDef(err,d,Vector{String},       :SNames                                     )
    d,err = checkType_setDef(err,d,Vector{String},       :XNames                                     )
    d,err = checkType_setDef(err,d,Int64,                :Nz                                         )
    d,err = checkType_setDef(err,d,Float64,              :V                                          )
    d,err = checkType_setDef(err,d,Float64,              :A                                          )
    d,err = checkType_setDef(err,d,Float64,              :LL,                 default=0.0            )
    d,err = checkType_setDef(err,d,Float64,              :Q,                                         ) 
    d,err = checkType_setDef(err,d,Vector{Float64},      :Xto                                        )
    d,err = checkType_setDef(err,d,Vector{Float64},      :Sto                                        )
    d,err = checkType_setDef(err,d,Vector{Float64},      :Pbo                                        )
    d,err = checkType_setDef(err,d,Vector{Float64},      :Sbo                                        )
    d,err = checkType_setDef(err,d,Float64,              :Lfo                                        )
    d,err = checkType_setDef(err,d,Array{Float64},       :Yxs                                        )
    d,err = checkType_setDef(err,d,Vector{Float64},      :Dt                                         )
    d,err = checkType_setDef(err,d,Vector{Float64},      :Db                                         )
    d,err = checkType_setDef(err,d,Vector{Float64},      :rho                                        )
    d,err = checkType_setDef(err,d,Float64,              :Kdet                                       ) 
    d,err = checkType_setDef(err,d,Tuple,                :Sin                                        )
    d,err = checkType_setDef(err,d,Tuple,                :srcS                                       )
    d,err = checkType_setDef(err,d,Tuple,                :mu                                         )
    d,err = checkType_setDef(err,d,Tuple,                :srcX                                       )
    d,err = checkType_setDef(err,d,Float64,              :Ptot,           default=sum(d[:Pbo])       )
    d,err = checkType_setDef(err,d,Int64,                :Nx,             default=length(d[:XNames]) )
    d,err = checkType_setDef(err,d,Int64,                :Ns,             default=length(d[:SNames]) )
        
    # Check if any errors 
    if err 
        println("")
        error("Issue with provided parameters
    (see https://markowkes.github.io/Biofilm.jl/parameters/) for more information.")
    end

    return d 
end

function checkType_setDef(err,d,type,name; default=nothing)
    myerr = false
    # Check if parameter exists or a default value is provided
    d,myerr = checkParamExists_setDef(myerr,d,name,default)
    # Check type
    if myerr == false
        try
            # Make Tuple of functions if needed 
            # (functions in Tuple will be tested later)
            if type == Tuple 
                if typeof(d[name]) <: Vector 
                    # convert to Tuple of functions
                    @reset d[name] = Tuple(d[name])
                elseif typeof(d[name]) <: Function
                    @reset d[name] = (d[name], )
                end
            end
            # Check type
            type(d[name])
        catch
            err = printError(err,"Parameter $name should be of type $type")
        end
    else
        err = true
    end
    return d,err
end

function checkParamExists_setDef(err,d,name,default)
    # Check if key missing from d
    if !haskey(d,name)
        # Check if default provided
        if default !== nothing 
            # Use default value 
            # @printf("Using default value :  %-19s = %s \n",name,default)
            @reset d[name] = default
        else
            err = printError(err,"Parameter $name missing from inputs!")
        end
    end
    return d,err
end

function printError(err,errMsg)
    try 
        error(errMsg)
    catch e
        printstyled(stderr,"ERROR: ", bold=true, color=:red)
        printstyled(stderr,sprint(showerror,e), color=:light_red)
        println(stderr)
    end
    err = true
    return err
end


function paramError(msg::Vararg{Any,N}) where {N}
    error("ISSUE WITH PARAMETERS:\n",
    "---------------------------------------------------------\n",
    "---------------------------------------------------------\n",
    string(msg...),
    "\n---------------------------------------------------------",
    "\n---------------------------------------------------------")
end

function checkParameters(p)

    @unpack Nx,Ns,Nz,Xto,Sto,Pbo,Sbo,Lfo,SNames,XNames,Title,mu,srcX,srcS,Sin,tFinal,outPeriod,plotPeriod,Yxs = p

    # Check provided initial conditions 
    Nx == length(Xto) || paramError("Number of Xto initial conditions should be Nx=", Nx)
    Ns == length(Sto) || paramError("Number of Sto initial conditions should be Ns=", Ns," found $(length(Sto))")
    Nx == length(Pbo) || paramError("Number of Pbo initial conditions should be Nx=", Nx)
    Ns == length(Sbo) || paramError("Number of Sbo initial conditions should be Ns=", Ns)

    # Check Names 
    Title isa String || paramError("Title should be a string describing the simulation.")
    SNames isa Array{String} && length(SNames) == Ns || paramError("SNames should be an array of Ns=",Ns," strings.")
    XNames isa Array{String} && length(XNames) == Nx || paramError("XNames should be an array of Nx=",Nx," strings.")

    # Check time parameters
    tFinal > 0 || paramError("Time should be a real number that is greater that zero")

    # Growthrate - check that can call it correctly
    for i in 1:Nx
        St=rand(Ns)
        Xt=rand(Nx)
        zm=0.5*Lfo
        try mu[i](St,Xt,Lfo,0.0,range(0,Lfo,length=Nz),p)
            muSize1 = size(mu[i](St,Xt,Lfo,0.0,zm,p),1)
            muSize2 = size(mu[i](St,Xt,Lfo,0.0,zm,p),2)
            muSize1 == 1 && muSize2 == 1 || 
            paramError("Error calling mu[",i,"].  mu returns an array of size (",muSize1,",",musize2," it should return an array of size (1,1).
                Check to make sure solutes are indexed correctly, e.g., S[1].")
            
        catch e
            paramError("Error calling mu[",i,"]. mu should be an array of Nx=",Nx," functions providing the growthrate of each particulate. 
            The inputs to each function should be (S,X,Lf,t,z,p) \n
                For example, if there are two particulates you might use:
                    mu=[(S,X,Lf,t,z,p) -> mumax * S[1,:] ./ (KM .+ S[1,:]), 
                        (S,X,Lf,t,z,p) -> mumax * S[2,:] ./ KM ],")
            println(e)
        end
    end

    # Source
    t=0.0
    for i in 1:Nx
        try srcX[i](Sto,Xto,Lfo,t,Lfo,p)
        catch e
            paramError("srcX should be an array of Nx=",Nx," functions providing the source of each particulate. 
            The inputs to each function should be (St,Xt,p) \n
                For example, if there are two particulates you might use:
                    srcX=[(S,X,Lf,t,z,p) -> -b*X[1],
                          (S,X,Lf,t,z,p) ->  b*X[1]],")
            println(e)
        end
    end

    for i in 1:Ns
        try srcS[i](Sto,Xto,Lfo,t,Lfo,p)
        catch e
            paramError("srcS should be an array of Ns=",Ns," functions providing the source of each solute. 
            The inputs to each function should be (St,Xt,t,z,p) \n
                For example, if there are two particulates you might use:
                    srcS=[(S,X,Lf,t,z,p) -> -b*X[1,:],
                          (S,X,Lf,t,z,p) ->  b*X[1,:]],")
            println(e)    
        end
    end

    # Sin
    for i in 1:Ns
        try Sin[i](0.0)
        catch e
            paramError("Sin should be an array of Ns=",Ns," functions providing the inflow solute concentrations. 
            The inputs to each function should be (t) \n
                For example, if there are two solutes you might use:
                    Sin=[(t) -> 25,
                         (t) -> 50],")
            println(e)
        end
    end

    # Yxs
    for i in 1:Nx
        try Yxs[i,1:Ns]
        catch e
            paramError("The size of Yxs is (",size(Yxs,1),",",size(Yxs,2),"). It should be (",Nx,",",Ns,") \n 
            For two solutes you might use:
                Yxs = [0.5 -0.7],  # Note the space between entries

            For two particulates you might use:
                Yxs = [0.5, -0.7], # Note the comma between entries
            or 
                Yxs= [ 0.5   # With a seperate line for each particulate
                      -0.7],
            ")
            println(e)
        end
    end

    # Replace any zero Yxs values wtih Inf to avoid division by zero
    Yxs[Yxs.==0.0].=Inf

    # Make sure plotPeriod is a multiple of outPeriod
    rem(plotPeriod,outPeriod) ≈ 0.0 || paramError("plotPeriod should be a multiple of outPerid.")
    
    return nothing
end