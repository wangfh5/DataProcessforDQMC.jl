using DataProcessforDQMC
using Test

@testset "DataProcessforDQMC.jl" begin
    include("statistics_format_value_error_test.jl")
    include("statistics_remove_outliers_test.jl")
    include("statistics_covariance_test.jl")
    include("scalar_measurements_iqrfence_kw_test.jl")
    include("structure_factor_outlier_filter_test.jl")
    include("correlation_ratio_outlier_filter_test.jl")
end
