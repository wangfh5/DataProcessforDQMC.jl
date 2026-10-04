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

    # Rounding up that carries into the next decade keeps the decimal place of the
    # unrounded error (PDG convention), with value precision aligned
    val_str, err_str = format_value_error(0.7012047252570854, 0.009466208599346684; format=:decimal)
    @test val_str == "0.701"
    @test err_str == "0.010"
    @test format_value_error(1.336759, 0.000979; format=:decimal) == ("1.3368", "0.0010")
    @test format_value_error(0.358305, 0.009802, 1; format=:decimal) == ("0.358", "0.010")
    @test format_value_error(1.336759, 0.000979) == ("1.3368e+00", "0.0010e0")
    # Two significant digits: 0.0996 -> 0.10 carries, shown as 0.100
    @test format_value_error(0.81, 0.0996, 2; format=:decimal) == ("0.810", "0.100")
    # No carry: unchanged
    @test format_value_error(1.3377, 0.00081, 1; format=:decimal) == ("1.3377", "0.0009")

    # Order-of-magnitude quoting (error_sig_digits=0): the error rounds up to the
    # position above its leading digit, and value precision stays aligned
    @test format_value_error(5.423756758475438, 0.0069434274377905714, 0; format=:decimal) == ("5.42", "0.01")
    # Idempotent for an error already sitting on a power of ten
    @test format_value_error(5.423756758475438, 0.01, 0; format=:decimal) == ("5.42", "0.01")
    # Scientific format
    @test format_value_error(5.423756758475438, 0.01, 0) == ("5.42e+00", "0.01e0")
end
