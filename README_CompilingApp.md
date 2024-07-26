# Steps to compile app

>> cd "directory containing Biofilm"
>> julia -q --project
julia> using PackageCompiler
julia> create_app("Biofilm","BiofilmCompiled",force=true,include_transitive_dependencies=false)

# Run App
>> BiofilmCompiled/bin/Biofilm /path/to/case/file/Case1.jl