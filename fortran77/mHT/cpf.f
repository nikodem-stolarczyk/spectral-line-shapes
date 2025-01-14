c----------------------------------------------------------------------
      function cpf_accurate(x, y)
c----------------------------------------------------------------------
c Computes the complex probability function using a rational
c series with 42 terms. Assumes Im(z) > 0 or Im(z) = 0.
c Input Parameters:
c   x : Real part of input complex parameter
c   y : Imaginary part of input complex parameter
c Output:
c   Complex probability function
c----------------------------------------------------------------------
      include 'constants.inc'
      double precision x, y
      double precision weidemann_constant_42, fft_constant_terms(37)
      integer i, number_of_fft_terms
      double complex cpf_accurate, z, z_ratio, cpf_z

      parameter (number_of_fft_terms = 37)
      parameter (weidemann_constant_42 = 5.449631621480024d0)

      data fft_constant_terms / -3.129493160727961d-14,
     &  -1.188364999909099d-14,  1.951777029849348d-13,
     &   1.790586243645278d-13, -1.184560208678836d-12,
     &  -2.069163661083667d-12,  6.430136110306704d-12,
     &   2.063579921011804d-11, -2.392389527320517d-11,
     &  -1.799169607159564d-10, -6.353807951660892d-11,
     &   1.282896083944607d-9,   2.636162411919059d-9,
     &  -5.468780625369738d-9,  -3.294773119114329d-8,
     &  -2.752070035718561d-8,   2.206733163926054d-7,
     &   8.511689670641750d-7,   4.936972061734341d-7,
     &  -6.617492208403963d-6,  -2.914574364851397d-5,
     &  -4.816473680511106d-5,   1.044072210002090d-4,
     &   1.070131083157417d-3,   4.631075611097791d-3,
     &   1.480296368764821d-2,   3.922970169744468d-2,
     &   9.038744880336540d-2,   1.857036333535562d-1,
     &   3.455278077566057d-1,   5.882708203344523d-1,
     &   9.230959991941070d-1,   1.342044484596932d0,
     &   1.814714451499866d0,    2.288734169675538d0,
     &   2.697763665856064d0,    2.975931371735470d0 /

      z = dcmplx(-y, x)
      z_ratio = (weidemann_constant_42 + z)/(weidemann_constant_42 - z)
      cpf_z = dcmplx(0.0d0, 0.0d0)

      do i = 1, number_of_fft_terms
         cpf_z = cpf_z + fft_constant_terms(i) * z_ratio** 
     &             (number_of_fft_terms - i)
      end do
      cpf_accurate = 2.0d0 * cpf_z  
     &             / (weidemann_constant_42 - z)**2.0d0
     &             + inverse_sqrt_pi / (weidemann_constant_42 - z)
      end function cpf_accurate
c----------------------------------------------------------------------
      function cpf_fast(x, y)
c----------------------------------------------------------------------
c Computes the complex probability function using Humlicek's algorithm 
c and rational series with 24 terms.
c Input Parameters:
c   x : Real part of input complex parameter
c   y : Imaginary part of input complex parameter
c Output:
c   Complex probability function
c----------------------------------------------------------------------
      include 'constants.inc'
      double precision x, y
      double precision hum1_threshold, weidemann_constant_24,
     &           fft_constant_terms(24)
      integer i, number_of_fft_terms
      double complex cpf_fast, z, z_ratio, cpf_z, t

      parameter (number_of_fft_terms = 24)
      parameter (hum1_threshold = 15.0d0)
      parameter (weidemann_constant_24 = 4.119534287814236d0)

      data fft_constant_terms / -1.513747622620502d-10,
     &   4.904820407381768d-9,   1.331045329581992d-9,
     &  -3.008282344381996d-8,  -1.912225887484805d-8,
     &   1.873834346505099d-7,   2.568264135399530d-7,
     &  -1.085647579417637d-6,  -3.038893184366094d-6,
     &   4.139461724429617d-6,   3.047106608295325d-5,
     &   2.433141546207148d-5,  -2.074843151143828d-4,
     &  -7.816642995626165d-4,  -4.936426901286291d-4,
     &   6.215006362949147d-3,   3.372336685531603d-2,
     &   1.083872348456673d-1,   2.654963959880772d-1,
     &   5.361139535729116d-1,   9.257087138588670d-1,
     &   1.394819673379119d0,    1.856286499205540d0,
     &   2.197858936531542d0 /

      if (abs(x) + y > hum1_threshold) then
         t = dcmplx(y, -x)
         cpf_z = inverse_sqrt_pi * t / (0.5d0 + t**2.0d0)
      else
         z = dcmplx(-y, x)
         z_ratio = (weidemann_constant_24 + z)
     &             / (weidemann_constant_24 - z)
         cpf_z = dcmplx(0.0d0, 0.0d0)

         do i = 1, number_of_fft_terms
            cpf_z = cpf_z + fft_constant_terms(i) * z_ratio**
     &                (number_of_fft_terms - i)
         end do

         cpf_z = 2.0d0 * cpf_z
     &         / (weidemann_constant_24 - z)**2.0d0
     &         + inverse_sqrt_pi / (weidemann_constant_24 - z)
      endif

      cpf_fast = cpf_z
      end function cpf_fast
