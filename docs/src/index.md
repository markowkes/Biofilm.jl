This package models the dynamics of a biofilm using the Julia programming language.  

# Installation 
```@contents
Pages = ["index.md"]
```
---

## Download & Install Julia

1. Download Julia from [julialang.org](https://julialang.org/downloads/)
2. Install Julia following the provided directions on the help pages
3. Launch Julia and you should see a terminal like the following (version number may be different)
![Julia REPL](images/juliaREPL.png)
4. If you see the terminal, you have successfully installed Julia!

!!! tip

    Try running a simple code like `julia> 5+5` to see that you can use Julia


## Add the Biofilm.jl package

- Run the following command in the Julia window

  ```
  using Pkg; Pkg.add(url="https://github.com/markowkes/Biofilm.jl")
  ```
  
!!! note

    This will download Biofilm.jl and all the dependencies and can take several minutes to finish

!!! tip

    If this command hangs on TriangularSolve, try using Julia 1.7.3


## Run Biofilm.jl 

Try running the examples.

1. Download ![Case1.jl](https://raw.githubusercontent.com/markowkes/Biofilm.jl/main/examples/Case1.jl) by right clicking and choosing *Save Link As*. Save the file to your Downloads folder (or other location of your choice)

2. Run the case by running the following within Julia on **Windows**

   ```
   using Biofilm   

   include("Downloads\\Case1.jl")  
   ```

   or within Julia on on **Mac/Linux**

   ```
   using Biofilm   

   include("Downloads/Case1.jl")  
   ```
    
!!! note

    The first time you run Biofilm.jl it will require time to compile. Subsequent runs should occur **much** faster. 

