function paramError(msg::Vararg{Any,N}) where {N}
    error("ISSUE WITH PARAMETERS:\n",
    "---------------------------------------------------------\n",
    "---------------------------------------------------------\n",
    string(msg...),
    "\n---------------------------------------------------------",
    "\n---------------------------------------------------------")
end

function checkParameters(p)

    @unpack Nx,Ns,Nz,Xto,Sto,Pbo,Sbo,Lfo,SNameSt,XtNames,Title,mu,srcX,srcS,Sin,tFinal,outPeriod,plotPeriod,Yxs = p

    # Check provided initial conditions 
    Nx == length(Xto) || paramError("Number of Xto initial conditions should be Nx=", Nx)
    Ns == length(Sto) || paramError("Number of Sto initial conditions should be Ns=", Ns)
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
        S=rand(Ns,Nz)
        X=rand(Nx,Nz)
        try mu[i](St,Xt,Lfo,0.0,range(0,Lfo,length=Nz),p)
            size(mu[i](St,Xt,Lfo,0.0,range(0,Lfo,length=Nz),p),1) == Nz &&
            size(mu[i](St,Xt,Lfo,0.0,range(0,Lfo,length=Nz),p),2) == 1  ||
            paramError("Error calling mu[",i,"].  mu returns an array of size ",muSize," it should return an array of size (1,",Nz,").
            Check to make sure substrates are indexed correctly, e.g., S[1,:].")
            
        catch e e
            paramError("Error calling mu[",i,"]. mu should be an array of Nx=",Nx," functions providing the growthrate of each particulate. 
            The inputs to each function should be (St,Xt,Lf,t,z,p) \n
                For example, if there are two particulates you might use:
                    mu=[(St,Xt,Lf,t,z,p) -> mumax * S[1,:] ./ (KM .+ S[1,:]), 
                        (St,Xt,Lf,t,z,p) -> mumax * S[2,:] ./ KM ],")
            println(e)
        end
    end

    # Source
    t=0.0
    for i in 1:Nx
        try srcX[i](Sto,Xto,t,p)
        catch e
            paramError("srcX should be an array of Nx=",Nx," functions providing the source of each particulate. 
            The inputs to each function should be (St,Xt,p) \n
                For example, if there are two particulates you might use:
                    srcX=[(S, X, t, p) -> -b*X[1,:],
                         (S, X, t, p) ->  b*X[1,:]],")
            println(e)
        end
    end

    for i in 1:Ns
        try srcS[i](Sto,Xto,t,p)
        catch e
            paramError("srcS should be an array of Ns=",Ns," functions providing the source of each substrate. 
            The inputs to each function should be (St,Xt,p) \n
                For example, if there are two particulates you might use:
                    srcS=[(S, X, t, p) -> -b*X[1,:],
                          (S, X, t, p) ->  b*X[1,:]],")
            println(e)    
        end
    end

    # Sin
    for i in 1:Ns
        try Sin[i](0.0)
        catch e
            paramError("Sin should be an array of Ns=",Ns," functions providing the inflow substrate concentrations. 
            The inputs to each function should be (t) \n
                For example, if there are two particulates you might use:
                    Sin=[(t) -> 25,
                         (t) -> 50],")
            println(e)
        end
    end

    # Yxs
    for i in 1:Nx
        try Yxs[i,1:Ns]
        catch e
            paramError("The size of Yxs is ",size(Yxs),". It should be (",Nx,",",Ns,")")
            println(e)
        end
    end

    # Replace any zero Yxs values wtih Inf to avoid division by zero
    Yxs[Yxs.==0.0].=Inf

    # Make sure plotPeriod is a multiple of outPeriod
    rem(plotPeriod,outPeriod) ≈ 0.0 || paramError("plotPeriod should be a multiple of outPerid.")
    
end