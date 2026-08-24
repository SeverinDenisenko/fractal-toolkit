program main_hurst
   use io, only: read_single_column
   use hurst, only: estimate_hurst_berg, estimate_hurst_yw, estimate_hurst_rs, estimate_hurst_lssd
   use version, only: max_ver_len, ver
   use precision, only: wp
   implicit none

   character(len=128) :: arg
   character(len=128) :: input
   character(len=max_ver_len) :: v
   real(wp), allocatable :: series(:)
   real(wp) :: H, a, H_err, a_err, sigma2
   integer :: i, m, ierr
   logical :: m_provided

   m_provided = .false.
   input = 'input.dat'

   i = 1
   do while (i <= command_argument_count())
      call get_command_argument(i, arg)

      select case (arg)
       case ('-v', '--version')
         call ver(v)
         print '(2a)', 'version ', v
         stop
       case ('-h', '--help')
         call print_help()
         stop
       case ('-i', '--input')
         if (i >= command_argument_count()) then
            print '(a)', 'Error: missing value for -i'
            stop 1
         end if
         i = i + 1
         call get_command_argument(i, arg)
         read(arg, *) input
       case ('-m', '--order')
         if (i >= command_argument_count()) then
            print '(a)', 'Error: missing value for -m'
            stop 1
         end if
         i = i + 1
         call get_command_argument(i, arg)
         read(arg, *) m
         m_provided = .true.
       case default
         print '(a,a,/)', 'Unrecognized command-line option: ', arg
         call print_help()
         stop 1
      end select
      i = i + 1
   end do

   call read_single_column(input, series)

   if (.not. m_provided) then
      m = int(sqrt(real(size(series), wp)))
   end if

   print '(a)', 'Given:'
   print '("     n = ", I8)', size(series)
   print '("     m = ", I8)', m

   call estimate_hurst_rs(series, H, H_err, sigma2)

   print '(a)', 'Rescaled Range analysis:'
   print '("     H = ", F6.3)', H
   print '("    dH = ", F6.3)', H_err
   print '("sigma2 = ", F6.3)', sigma2

   call estimate_hurst_berg(series, m, H, a, H_err, a_err, sigma2)

   print '(a)', 'Burg method:'
   print '("     H = ", F6.3)', H
   print '("     a = ", F6.3)', a
   print '("    dH = ", F6.3)', H_err
   print '("    da = ", F6.3)', a_err
   print '("sigma2 = ", F6.3)', sigma2

   call estimate_hurst_yw(series, m, H, a, H_err, a_err, sigma2)

   print '(a)', 'Yule-Walker method:'
   print '("     H = ", F6.3)', H
   print '("     a = ", F6.3)', a
   print '("    dH = ", F6.3)', H_err
   print '("    da = ", F6.3)', a_err
   print '("sigma2 = ", F6.3)', sigma2

   call estimate_hurst_lssd(series, 1, 50, H, H_err, ierr=ierr)

   if (ierr == 0) then
      print '(a)', 'Koutsoyiannis (LSSD) method:'
      print '("     H = ", F6.3)', H
      print '("    dH = ", F6.3)', H_err
   end if

contains
   subroutine print_help()
      print '(a)', 'Estimate Hurst exponent'
      print '(a)', ''
      print '(a)', 'usage: hurst [options]'
      print '(a)', ''
      print '(a)', 'cmdline options:'
      print '(a)', ''
      print '(a)', '  -v, --version     print version information and exit'
      print '(a)', '  -h, --help        print usage information and exit'
      print '(a)', '  -i, --input       select input file'
      print '(a)', '  -m, --order       select berg AR filter length'
   end subroutine print_help
end program main_hurst
