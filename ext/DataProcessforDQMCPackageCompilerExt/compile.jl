using DataProcessforDQMC: DataProcessforDQMC, Algorithm
using PackageCompiler: PackageCompiler

function DataProcessforDQMC.compile(
  ::Algorithm{:PackageCompiler};
  dir::AbstractString=DataProcessforDQMC.default_compile_dir(),
  filename::AbstractString=DataProcessforDQMC.default_compile_filename(),
)
  if !isdir(dir)
    println("""The directory "$dir" doesn't exist yet, creating it now.""")
    println()
    mkdir(dir)
  end
  path = joinpath(dir, filename)
  println(
    """Creating the system image "$path" containing the compiled version of DataProcessforDQMC. This may take a few minutes.""",
  )
  PackageCompiler.create_sysimage(
    :DataProcessforDQMC;
    sysimage_path=path,
    precompile_execution_file=joinpath(@__DIR__, "precompile_dataprocessfordqmc.jl"),
  )
  println(DataProcessforDQMC.compile_note(; dir, filename))
  return path
end