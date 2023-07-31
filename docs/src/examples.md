
# Examples
```@contents
Pages = ["examples.md"]
```
---

## Case 1 - Single Solute and Particulate Case
 **Download link: [Case1.jl](https://raw.githubusercontent.com/markowkes/Biofilm.jl/main/examples/Case1.jl)** (right-click & Save Link As)

This case is the same the case shown on the [Installation](@ref) page.  It simulates a single particulate ("Heterotroph") and a single solute ("Nutrient").  The particulate has a growthrate of 
```math
\mu = \mu_\mathrm{max} \frac{S}{K_M + S}
```

## Case 2 - Multiple Solutes
**Download link: [Case2.jl](https://raw.githubusercontent.com/markowkes/Biofilm.jl/main/examples/Case2.jl)** (right-click & Save Link As)

This is a simple example of how multiple solutes can be simulated.  Solute 1 is used by "Bug" while solute 2 is completely indpendent (not very interesting).

## Case 3 - Live/Dead Bugs 
**Download link: [Case3.jl](https://raw.githubusercontent.com/markowkes/Biofilm.jl/main/examples/Case3.jl)** (right-click & Save Link As)

This example has living bugs that use the solute to grow.  The bugs die and the concentration of dead bugs is also computed. The source term, **src**, is used to transfer living bugs to dead bugs.
```math
\mathrm{src}_{X_1} = \mathrm{src}_\mathrm{living} = -b X_1\\
\mathrm{src}_{X_2} = \mathrm{src}_\mathrm{dead} = b X_1
```
where `b=0.1` is a constant.  The first source term reduces living bugs and the second adds these to the dead bugs. 

## Case 4 - Multiple Particulate and Solutes
**Download link: [Case4.jl](https://raw.githubusercontent.com/markowkes/Biofilm.jl/main/examples/Case4.jl)** (right-click & Save Link As)

## Case 5 - Phototroph
**Download link: [Case5_Photroph.jl](https://raw.githubusercontent.com/markowkes/Biofilm.jl/main/examples/Case5_Photroph.jl)** (right-click & Save Link As)

This case shows how you can simulate a species that grows in light.  The light turns off and on throughout each day and also only penetrates the top of the biofilm.  Since the light has discontinuities, the case has an additional parameter
```
    discontinuityPeriod=0.25,
```
that tells the solver the discontinuities can occur every 0.25 days. 


