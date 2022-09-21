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

## Download Biofilm.jl

- Download the Biofilm.jl package using either of the methods below
    1. [Direct download of current version](https://github.com/markowkes/Biofilm.jl/archive/refs/heads/main.zip) - use this if you do not have experience with Git.
        - Unzip the package in a location you can find later
    2. [Git repository](https://github.com/markowkes/Biofilm.jl) - access repository including package updates
        - Clone repository into a location you can find later
!!! tip

    Put the Biofilm.jl package in an easy to remember location such as `user/Biofilm`

## Run Biofilm.jl
  - Start Julia 
  
  - Change directory to the location you have Biofilm.jl stored using the `cd("PATH")` command.  For example 

    ```
    julia> cd("user/Biofilm")
    ```
    Note that `julia>` is the Julia prompt and is not part of the command. 

  - Run one of the examples located in the examples directory.  For example, Case 1 can be run using 

    ```
    julia> include("examples/Case1.jl")
    ```
    
!!! note

    The first time you run Biofilm.jl it will require time to download and compile dependencies.  Subsequent runs should occur much faster. 

