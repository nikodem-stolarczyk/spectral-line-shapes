      program example_cpf
      implicit none
      include 'constants.inc'
      
      ! Declare variables
      double precision cpf_input_real, cpf_input_imag
      double complex cpf_accurate, cpf_fast
      double complex cpf_accurate_output, cpf_fast_output

      ! Initialize example parameters for the CPF functions
      cpf_input_real = 1.0d0  ! dimensionless
      cpf_input_imag = 1.0d0  ! dimensionless

      ! Call cpf_accurate
      cpf_accurate_output = cpf_accurate(cpf_input_real, cpf_input_imag)

      ! Call cpf_fast
      cpf_fast_output = cpf_fast(cpf_input_real, cpf_input_imag)

      ! Output results
      write(*, '(A)') 'the output of the cpf_accurate function: '
      write(*, '(2F22.15)') real(cpf_accurate_output),
     &                        imag(cpf_accurate_output)

      write(*, '(A)') 'the output of the cpf_fast function:     '
      write(*, '(2F22.15)') real(cpf_fast_output), imag(cpf_fast_output)

      end program example_cpf
