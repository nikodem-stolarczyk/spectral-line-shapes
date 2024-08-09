program example_cpf

   use, intrinsic :: iso_fortran_env, only: int32, dp => real64
   use cpf_module, only: cpf_accurate, cpf_fast
   implicit none
   ! example parameters for the cpf functions
   real(dp)    :: cpf_input_real = 1.0_dp ! dimensionless
   real(dp)    :: cpf_input_imag = 1.0_dp ! dimensionless
   complex(dp) :: cpf_accurate_output ! Example output of cpf_accurate(cpf_input_real, cpf_input_imag)
   complex(dp) :: cpf_fast_output ! Example output of cpf_fast(cpf_input_real, cpf_input_imag)
   cpf_accurate_output = cpf_accurate(cpf_input_real, cpf_input_imag)
   cpf_fast_output = cpf_fast(cpf_input_real, cpf_input_imag)
   
   write(*, '(A, 2F22.15)') 'the output of the cpf_accurate function: ',&
      real(cpf_accurate_output), aimag(cpf_accurate_output)
   write(*, '(A, 2F22.15)') 'the output of the cpf_fast function:     ',&
      real(cpf_fast_output), aimag(cpf_fast_output)

end program example_cpf
