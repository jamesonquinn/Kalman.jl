using Pkg
Pkg.add("LinearAlgebra")
Pkg.add("Distributions")
Pkg.add("StatsBase")
Pkg.add("Memoize")
Pkg.add("Compat")
Pkg.add("ForwardDiff")
Pkg.add("DataFrames")
Pkg.add("CSV")
Pkg.add("DataStructures")
Pkg.add("Revise")
#Pkg.add("Plots")
#Pkg.add("Plotly")
#Pkg.add("Gadfly")
Pkg.add("NaNMath")
Pkg.add("SharedArrays")

if false
  using Distributed
  @sync @distributed for i in 1:10
    print(i)
    print("\n")
  end

  print("\n\n\n")
end
