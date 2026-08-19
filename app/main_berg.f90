program main_berg
   use stdlib_error, only: check
   use io, only: read_single_column
   use hurst, only: estimate_hurst_berg
   use version, only: max_ver_len, ver
   use precision, only: wp
   implicit none

   character(len=128) :: arg
   character(len=128) :: input
   character(len=max_ver_len) :: v
   real(wp), allocatable :: series(:)
   real(wp) :: H, a, H_err, a_err, sigma2
   integer :: i, m

   input = 'input.dat'
   m = 3

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
       case ('-m', '--length')
         if (i >= command_argument_count()) then
            print '(a)', 'Error: missing value for -m'
            stop 1
         end if
         i = i + 1
         call get_command_argument(i, arg)
         read(arg, *) m
       case default
         print '(a,a,/)', 'Unrecognized command-line option: ', arg
         call print_help()
         stop 1
      end select
      i = i + 1
   end do

   call read_single_column(input, series)

   call estimate_hurst_berg(series, m, H, a, H_err, a_err, sigma2)

   print '("H =      ", F6.3)', H
   print '("a =      ", F6.3)', a
   print '("dH =     ", F6.3)', H_err
   print '("da =     ", F6.3)', a_err
   print '("sigma2 = ", F6.3)', sigma2
contains
   subroutine print_help()
      print '(a)', 'usage: berg [options]'
      print '(a)', ''
      print '(a)', 'cmdline options:'
      print '(a)', ''
      print '(a)', '  -v, --version     print version information and exit'
      print '(a)', '  -h, --help        print usage information and exit'
      print '(a)', '  -i, --input       select input file'
      print '(a)', '  -m, --length      select berg AR filter length'
   end subroutine print_help
end program main_berg
