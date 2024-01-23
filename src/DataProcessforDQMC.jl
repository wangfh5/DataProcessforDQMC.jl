module DataProcessforDQMC

println("Welcome to DataProcessforDQMC!")

export myfunction1, myfunction2 # export functions
# Write your package code here.
myfunction1(x) = x + 1
myfunction2(x) = x * 2

include("format-dataframes.jl")

end
