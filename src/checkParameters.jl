function paramError(msg::Vararg{Any,N}) where {N}
    error("ISSUE WITH PARAMETERS:\n",
    "---------------------------------------------------------\n",
    "---------------------------------------------------------\n",
    string(msg...),
    "\n---------------------------------------------------------",
    "\n---------------------------------------------------------")
end

function checkParameters(p)

    @unpack Nx,Ns,Nz,Xo,So,Pbo,Sbo,Lfo,SNames,XNames,Title,mu,src,Sin,tFinal,outPeriod,Yxs = p

    # Check provided initial conditions 
    Nx == length(Xo)  || paramError("Number of Xo initial conditions should be Nx=", Nx)
    Ns == length(So)  || paramError("Number of So initial conditions should be Ns=", Ns)
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
        try mu[i](So,Xo,Lfo,0.0,range(0,Lfo,length=Nz),p)
        catch
            paramError("mu should be an array of Nx=",Nx," functions providing the growthrate of each particulate. 
            The inputs to each function should be (S,X,Lf,t,z,p) \n
                For example, if there are two particulates you might use:
                    mu=[(S,X,Lf,t,z,p) -> mumax * S[1,:] ./ (KM .+ S[1,:]), 
                        (S,X,Lf,t,z,p) -> mumax * S[2,:] ./ KM ],")
        end
    end

    # Source
    for i in 1:Nx
        try src[i](So,Xo,p)
        catch
            paramError("src should be an array of Nx=",Nx," functions providing the source of each particulate. 
            The inputs to each function should be (S,X,p) \n
                For example, if there are two particulates you might use:
                    src=[(S, X, p) -> -b*X[1,:],
                         (S, X, p) ->  b*X[1,:]],")
        end
    end

    # Sin
    for i in 1:Ns
        try Sin[i](0.0)
        catch
            paramError("Sin should be an array of Ns=",Ns," functions providing the inflow substrate concentrations. 
            The inputs to each function should be (t) \n
                For example, if there are two particulates you might use:
                    Sin=[(t) -> 25,
                         (t) -> 50],")
        end
    end

    # Yxs
    for i in 1:Nx
        try Yxs[i,1:Ns]
        catch
            paramError("The size of Yxs is ",size(Yxs),". It should be (",Nx,",",Ns,")")
        end
    end

    # Replace any zero Yxs values wtih Inf to avoid division by zero
    Yxs[Yxs.==0.0].=Inf
end