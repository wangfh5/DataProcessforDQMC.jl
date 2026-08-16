using Test
using DataProcessforDQMC: format_value_error

@testset "format_value_error examples" begin
    @test format_value_error(2.36738, 0.0023) == ("2.367e+00", "0.003e0")
    @test format_value_error(2367.38, 23, 2) == ("2.367e+03", "0.023e3")
    @test format_value_error(2.36738, 0.0023; format=:decimal) == ("2.367", "0.003")
    @test format_value_error(2367.38, 23; format=:decimal) == ("2370", "30")

    # Regression: error==0 should NOT round value to integer (e.g. 0.5 -> 0)
    @test format_value_error(0.5, 0.0) == ("5e-01", "0e-1")
    @test format_value_error(0.5, 0.0; format=:decimal) == ("0.5", "0")

    # Regression: rounding up a small error should keep value precision aligned
    val_str, err_str = format_value_error(0.7012047252570854, 0.009466208599346684; format=:decimal)
    @test val_str == "0.70"
    @test err_str == "0.01"

    # Order-of-magnitude quoting (error_sig_digits=0): the error rounds up to the
    # position above its leading digit, and value precision stays aligned
    @test format_value_error(5.423756758475438, 0.0069434274377905714, 0; format=:decimal) == ("5.42", "0.01")
    # Idempotent for an error already sitting on a power of ten
    @test format_value_error(5.423756758475438, 0.01, 0; format=:decimal) == ("5.42", "0.01")
    # Scientific format
    @test format_value_error(5.423756758475438, 0.01, 0) == ("5.42e+00", "0.01e0")
end
