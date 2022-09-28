This package models the dynamics of a biofilm using the Julia programming language.  

# Installation 
```@contents
Pages = ["index.md"]
```
---

## Download & Install Julia

1. Download the latest version of Julia from [julialang.org](https://julialang.org/downloads/)

2. Install Julia following the provided directions on the [help page](https://julialang.org/downloads/platform/).  The default installation directory should be fine.  On Windows be sure to select the option to add Julia to the path to allow VS Code (next section) to find Julia.  On Mac or Linux be sure to follow the instructions on adding Julia to PATH.

3. Launch Julia by finding the icon where you installed the program.  If you have been successful, you should see the REPL like this image (version number may be different)
![Julia REPL](images/juliaREPL.png)

!!! tip

    Try running a simple code like `julia> 5+5` to see that you can use Julia.  You could even try learning more about the language by doing a tutorial like this one: "[From zero to Julia!](https://techytok.com/from-zero-to-julia/)". 

## VS Code - GUI for Julia
Julia is a programming languages and you can use it in many ways.  One popular way to run Julia code is through Visual Studio Code, which allows you to edit, run, and see results from a code.

1. Download, install, and open VS Code.  You can find the installers on [code.visualstudio.com](https://code.visualstudio.com/download)
   
2. Install the Julia extension in VS Code by selecting `File [Code on Mac]` >> `Preferences` >> `Extensions`, then search for Julia.  You should see an extension called Julia provided by julialang.  Click the `Install` button.
   ![Julia Extension](images/juliaExtension.png)

3. Start the Julia REPL in VS Code by opening the command-pallet using `Ctrl-Shift-p` on Windows or `Cmd-Shfit-p` on Mac then searching for and running `Julia: Start REPL`.  This will start Julia inside of VS Code.


## Add the Biofilm.jl package

- Run the following command in the Julia REPL within VS Code.

  ```
  using Pkg; Pkg.add(url="https://github.com/markowkes/Biofilm.jl")
  ```
This will download Biofilm.jl and all the dependencies and can take several minutes to finish.


## Run Biofilm.jl 

Try running the examples.

1. Download [Case1.jl](https://raw.githubusercontent.com/markowkes/Biofilm.jl/main/examples/Case1.jl) by right clicking and choosing *Save Link As*. Save the file to your Downloads folder (or other location of your choice)

2. Run the case by opening the file in VS Code and clicking the play button in the top right corner.  You should see output to the REPL and a plot like this:
   ![Case 1 Run Button](images/Case1_Run.png)
The top row shows the biomass (particulate) and substrate concentrations and the biofilm thickness as a function of time.  The bottom row shows the particulate volume fraction and substrate concentrations as a function of position within the biofilm at the output time (t=1.0 as indicated by the title).
   
3. Try editing the case file or explore other [Examples](@ref)
    
!!! note

    The first time you run Biofilm.jl it will require time to compile. Subsequent runs should occur **much** faster. 
