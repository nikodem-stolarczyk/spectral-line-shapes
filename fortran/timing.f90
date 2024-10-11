program cpf_test
   use, intrinsic :: iso_fortran_env, only: int32, dp => real64
   use cpf_module, only: cpf_accurate, cpf_fast
!----------------------------------------------------------------------!
   implicit none
!----------------------------------------------------------------------!   
   complex(dp) :: cpf_accurate_output
      ! Example output of cpf_accurate(cpf_input_real, cpf_input_imag)
   complex(dp) :: cpf_fast_output
      ! Example output of cpf_fast(cpf_input_real, cpf_input_imag)
   character(:), allocatable :: filename_r, filename_i
   integer(int32), parameter :: grid = 1000
   real(dp), dimension(grid) :: xprep, x, y
   real(dp), dimension(grid*grid) :: weid32r, weid32i
!----------------------------------------------------------------------!
   integer(int32) :: ii, ix, iy
   real(dp) :: t1, t2, t_total, step, start_, end_
!----------------------------------------------------------------------!
   filename_r = 'weid24r_f.csv'
   filename_i = 'weid24r_i.csv'
!----------------------------------------------------------------------!
   print '(a,1x,i0,a)',                                               &
      'Evaluating complex error function on the', grid, '-point grid...'
   print '(a,1x,a,1x,a,1x,a,1x,a)',                                    &
            'The results will be written to', filename_r, 'and',       &
            filename_i, 'files'
!----------------------------------------------------------------------!
   start_ = 0
   end_ = 3
! Generate xprep array using linspace logic
  do ii = 1, grid
     xprep(ii) = -8.0d0 + (2.0d0 - (-8.0d0)) * (ii - 1) / (grid - 1)
  end do
   x = 10.0_dp ** xprep
  ! Generate y array using linspace logic
  do ii = 1, grid
     y(ii) = start_ + (end_ - start_) * (ii - 1) / (grid - 1)
  end do
!----------------------------------------------------------------------!
   call cpu_time(t1)
!----------------------------------------------------------------------!
   ii = 1
   do ix = 1, grid
      do iy = 1, grid
         cpf_accurate_output = cpf_accurate(y(iy), x(ix))
         ii = ii + 1
      enddo
   enddo
   call cpu_time(t2)
!----------------------------------------------------------------------!
   t_total = t2 - t1
   print '(a,1x,f10.4,1x,a)', 'cpf_accurate time', t_total, 's'
!----------------------------------------------------------------------!
   call cpu_time(t1)
!----------------------------------------------------------------------!
  ii = 1
   do ix = 1, grid
      do iy = 1, grid
         cpf_fast_output = cpf_fast(y(iy), x(ix))
         ii = ii + 1
      enddo
   enddo
   call cpu_time(t2)
!----------------------------------------------------------------------!
   t_total = t2 - t1
   print '(a,1x,f10.4,1x,a)', 'cpf_fast time', t_total, 's'
!!----------------------------------------------------------------------!
!   open(unit = 11, file = trim(filename_r), action = 'write')
!   open(unit = 12, file = trim(filename_i), action = 'write')
!!----------------------------------------------------------------------!
!! caution: unformatted output, keeps 16 significant digits
!!----------------------------------------------------------------------!
!   do ii = 1, size(weid32r)
!      write(11, *) weid32r(ii)
!      write(12, *) weid32i(ii)
!   enddo
!   close(11)
!   close(12)
!----------------------------------------------------------------------!
end program cpf_test
