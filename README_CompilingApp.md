# Steps to compile app
```
>> cd "directory containing Biofilm"
>> julia -q --project
julia> using PackageCompiler
julia> create_app("Biofilm","BiofilmCompiled",force=true,include_transitive_dependencies=false)
```

# Run App using a Case file 
- any in the test directory should work
- the ones in the examples directory also work but call the solver twice
For example
```
>> ./BiofilmCompiled/bin/Biofilm Biofilm/test/Case1.jl
```
