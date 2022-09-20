
# Examples
```@contents
Pages = ["examples.md"]
```
---

## Case 1 - Single Substrate and Particulate Case

### Running the case
  1. Start Julia and change the directory to location you installed Julia

  2. Run the case using 

     ```
     include("examples/Case1.jl")
     ```

  3. The simulation runs and should produce an output like ![Case 1 Output](images/Case1_final.svg)  The top row shows the biomass (particulate) and substrate concentrations and the biofilm thickness as a function of time.  The bottom row shows the particulate volume fraction and substrate concentrations as a function of position within the biofilm at the output time (t=1.0 as indicated by the title).