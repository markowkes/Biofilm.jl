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
  julia> using Pkg; Pkg.add(url="https://github.com/markowkes/Biofilm.jl")
  ```
  Note that `julia>` is the Julia prompt and is not part of the command.  This will download Biofilm.jl and all the dependencies


## Run Biofilm.jl

  - Run one of the examples located in the examples directory.  For example, Case 1 can be run using 

    ```
    julia> include("examples/Case1.jl")
    ```
    
!!! note

    The first time you run Biofilm.jl it will require time to download and compile dependencies.  Subsequent runs should occur much faster. 

