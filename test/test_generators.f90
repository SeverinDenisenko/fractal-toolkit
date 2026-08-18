program generte_test
   use precision, only: wp
   use generators, only: generte_fgn, generte_white, generte_fgn_integrate
   use stdlib_error, only: check
   implicit none

   real(wp) :: result(5)

   call generte_fgn(result)
   call check(maxval(abs(result - [-0.65635862872429629, 1.6005188979418756, -1.0128788452223753E-002, 1.3881885355976282, 1.1373458198637005])) < 1e-5_wp)

   call generte_white(result)
   call check(maxval(abs(result - [-0.1911449112858024_wp, -0.2687183825051316_wp, -0.0014389361844231_wp, -0.3206724268956433_wp, -0.3209329317255407_wp])) < 1e-5_wp)

   call generte_fgn_integrate(result, 1.0_wp)
   call check(maxval(abs(result - [-0.65635862872429629, 0.94416026921757934, .93403148076535558, 2.3222200163629836, 3.4595658362266839])) < 1e-5_wp)
end program generte_test
