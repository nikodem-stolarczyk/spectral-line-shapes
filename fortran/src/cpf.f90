module cpf_module
!----------------------------------------------------------------------!
! This module provides 2 functions that compute complex probability
! function (CPF):
! 1. cpf_accurate - computes the CPF using rational series with 42 terms
! 2. cpf_fast - computes the CPF using Humlicek's algorithm and rational
!    series with 24 terms
!----------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, dp => real64
   implicit none
!----------------------------------------------------------------------!
   real(dp), parameter :: inverse_sqrt_pi = 0.5641895835477563_dp
!----------------------------------------------------------------------!
   contains
!----------------------------------------------------------------------!
   function cpf_accurate(x, y) result(cpf_z)
      !----------------------------------------------------------------!
      ! Computes the complex probability function using a rational
      ! series with 42 terms. It is assumed that Im(z) > 0 or Im(z) = 0.
      ! A series was simplified to 37 terms introducing less than
      ! 10^(-17) deviations on mHT profile.
      !
      ! Input/Output Parameters of Routine
      ! ---------------------------------------------------------------!
      ! x : Real part of input complex parameter
      ! y : Imaginary part of input complex parameter
      !
      ! The function provides one output:
      ! ---------------------------------------------------------------!
      ! (1): Complex probability function
      !----------------------------------------------------------------!
      real(dp), intent(in) :: x, y
      complex(dp) :: cpf_z
      !----------------------------------------------------------------!
      integer(int32), parameter :: number_of_fft_terms = 37
      real(dp), parameter :: weidemann_constant_42 = 5.449631621480024_dp
      real(dp), dimension(number_of_fft_terms) :: fft_constant_terms=(/&
         -3.129493160727961e-14_dp, -1.188364999909099e-14_dp,         &
          1.951777029849348e-13_dp,  1.790586243645278e-13_dp,         &
         -1.184560208678836e-12_dp, -2.069163661083667e-12_dp,         &
          6.430136110306704e-12_dp,  2.063579921011804e-11_dp,         &
         -2.392389527320517e-11_dp, -1.799169607159564e-10_dp,         &
         -6.353807951660892e-11_dp,  1.282896083944607e-9_dp,          &
          2.636162411919059e-09_dp, -5.468780625369738e-9_dp,          &
         -3.294773119114329e-8_dp,  -2.752070035718561e-8_dp,          &
          2.206733163926054e-7_dp,   8.511689670641750e-7_dp,          &
          4.936972061734341e-7_dp,  -6.617492208403963e-6_dp,          &
         -2.914574364851397e-5_dp,  -4.816473680511106e-5_dp,          &
          1.044072210002090e-4_dp,   1.070131083157417e-3_dp,          &
          4.631075611097791e-3_dp,   1.480296368764821e-2_dp,          &
          3.922970169744468e-2_dp,   9.038744880336540e-2_dp,          &
          1.857036333535562e-1_dp,   3.455278077566057e-1_dp,          &
          5.882708203344523e-1_dp,   9.230959991941070e-1_dp,          &
          1.342044484596932_dp,      1.814714451499866_dp,             &
          2.288734169675538_dp,      2.697763665856064_dp,             &
          2.975931371735470_dp/)
      !----------------------------------------------------------------!
      complex(dp) :: z, z_ratio
      integer(int32) :: i
      !----------------------------------------------------------------!
      z       = cmplx(-y, x, kind=dp)
      z_ratio = (weidemann_constant_42 + z) / (weidemann_constant_42 - z)
      cpf_z   = (0.0_dp, 0.0_dp)
      do i = 1, number_of_fft_terms
         cpf_z = cpf_z + fft_constant_terms(i)                         &
            * z_ratio**(real(number_of_fft_terms - i, kind = dp))
      end do
      cpf_z = 2.0_dp * cpf_z / (weidemann_constant_42 - z)**2.0_dp     &
         + inverse_sqrt_pi / (weidemann_constant_42 - z)
      !----------------------------------------------------------------!
   end function cpf_accurate
!----------------------------------------------------------------------!
   function cpf_fast(x, y) result(cpf_z)
      !----------------------------------------------------------------!
      ! Computes the complex probability function using Humlicek's 
      ! algorithm in its first subregion (Source: 10.1016/0022-4073(82)90078-4) 
      ! and using a rational series with 24 terms in other subregions.
      !
      ! Input/Output Parameters of Routine
      ! ---------------------------------------------------------------!
      ! x : Real part of input complex parameter
      ! y : Imaginary part of input complex parameter
      !
      ! The function has one output:
      ! ---------------------------------------------------------------!
      ! (1): Complex probability function
      ! ---------------------------------------------------------------!
      real(dp), intent(in) :: x, y
      complex(dp) :: cpf_z
      ! ---------------------------------------------------------------!
      integer(int32), parameter :: number_of_fft_terms = 24
      real(dp), parameter :: hum1_threshold = 15.0_dp
      real(dp), parameter :: weidemann_constant_24 = 4.119534287814236_dp
      real(dp), dimension(number_of_fft_terms) :: fft_constant_terms=(/&
         -1.513747622620502e-10_dp, 4.904820407381768e-9_dp,           &
          1.331045329581992e-9_dp, -3.008282344381996e-8_dp,           &
         -1.912225887484805e-8_dp,  1.873834346505099e-7_dp,           &
          2.568264135399530e-7_dp, -1.085647579417637e-6_dp,           &
         -3.038893184366094e-6_dp,  4.139461724429617e-6_dp,           &
          3.047106608295325e-5_dp,  2.433141546207148e-5_dp,           &
         -2.074843151143828e-4_dp, -7.816642995626165e-4_dp,           &
         -4.936426901286291e-4_dp,  6.215006362949147e-3_dp,           &
          3.372336685531603e-2_dp,  1.083872348456673e-1_dp,           &
          2.654963959880772e-1_dp,  5.361139535729116e-1_dp,           &
          9.257087138588670e-1_dp,  1.394819673379119_dp,              &
          1.856286499205540_dp,     2.197858936531542_dp/)
      !----------------------------------------------------------------!
      complex(dp) :: t, z, z_ratio
      integer(int32) :: i
      !----------------------------------------------------------------!
      if (abs(x) + y > hum1_threshold) then
         !-------------------------------------------------------------!
         t = cmplx(y, -x, kind=dp)
         cpf_z = inverse_sqrt_pi * t / (0.5_dp + t**2.0_dp)
         !-------------------------------------------------------------!
      else
         !-------------------------------------------------------------!
         z       = cmplx(-y, x, kind=dp)
         z_ratio = (weidemann_constant_24 + z) / (weidemann_constant_24 - z)
         cpf_z   = (0.0_dp, 0.0_dp)
         do i = 1, number_of_fft_terms
            cpf_z = cpf_z + fft_constant_terms(i)                      &
               * Z**(real(number_of_fft_terms - i, kind = dp))
         end do
         cpf_z = 2.0_dp * cpf_z / (weidemann_constant_24 - z)**2.0_dp  &
            + inverse_sqrt_pi / (weidemann_constant_24 - z)
         !-------------------------------------------------------------!
      endif
      !----------------------------------------------------------------!
  end function cpf_fast
!----------------------------------------------------------------------!
end module cpf_module
