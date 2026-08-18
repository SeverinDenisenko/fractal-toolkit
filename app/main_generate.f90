program main_generate
   use stdlib_error, only: check
   use generators, only: generte_white_integrate, generte_fgn_integrate, default_seed
   use io, only: write_table_file
   use precision, only: wp
   use version, only: max_ver_len, ver
   implicit none

   character(len=128) :: arg
   character(len=128) :: output
   character(len=max_ver_len) :: v
   real(wp), allocatable :: series(:)
   integer :: n, i, s
   real(wp) :: a

   output = 'output.dat'
   a = 0.5_wp
   n = 10
   s = default_seed

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
       case ('-n', '--length')
         if (i >= command_argument_count()) then
            print '(a)', 'Error: missing value for -n'
            stop 1
         end if
         i = i + 1
         call get_command_argument(i, arg)
         read(arg, *) n
       case ('-a', '--intorder')
         if (i >= command_argument_count()) then
            print '(a)', 'Error: missing value for -a'
            stop 1
         end if
         i = i + 1
         call get_command_argument(i, arg)
         read(arg, *) a
       case ('-s', '--seed')
         if (i >= command_argument_count()) then
            print '(a)', 'Error: missing value for -s'
            stop 1
         end if
         i = i + 1
         call get_command_argument(i, arg)
       case ('-o', '--output')
         if (i >= command_argument_count()) then
            print '(a)', 'Error: missing value for -o'
            stop 1
         end if
         i = i + 1
         call get_command_argument(i, arg)
         read(arg, *) output
       case default
         print '(a,a,/)', 'Unrecognized command-line option: ', arg
         call print_help()
         stop 1
      end select
      i = i + 1
   end do

   allocate(series(n))
   call generte_fgn_integrate(series, a, s)
   call write_table_file(output, series)
contains
   subroutine print_help()
      print '(a)', 'usage: generate [options]'
      print '(a)', ''
      print '(a)', 'cmdline options:'
      print '(a)', ''
      print '(a)', '  -v, --version     print version information and exit'
      print '(a)', '  -h, --help        print usage information and exit'
      print '(a)', '  -n, --length      select series length'
      print '(a)', '  -a, --intorder    select integration order'
      print '(a)', '  -s, --seed        select seed'
      print '(a)', '  -o, --output      set output'
   end subroutine print_help
end program main_generate
