# Steps to compile app

>> cd "directory containing Biofilm"
>> julia -q --project
julia> using PackageCompiler
julia> create_app("Biofilm","BiofilmCompiled",force=true,include_transitive_dependencies=false)
