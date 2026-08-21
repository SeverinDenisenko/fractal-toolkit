program main_fracdiff
   use io, only: read_single_column, write_table_file
   use frac, only: frac_diff_simple
   use version, only: max_ver_len, ver
   use precision, only: wp
   implicit none

   character(len=128) :: arg
   character(len=128) :: input, output
   character(len=max_ver_len) :: v
   real(wp), allocatable :: series(:), result(:)
   integer :: i
   real(wp) :: a

   input = 'input.dat'
   output = 'output.dat'
   a = 1

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
       case ('-o', '--output')
          if (i >= command_argument_count()) then
             print '(a)', 'Error: missing value for -o'
            stop 1
         end if
         i = i + 1
         call get_command_argument(i, arg)
         read(arg, *) output
       case ('-a', '--order')
         if (i >= command_argument_count()) then
            print '(a)', 'Error: missing value for -a'
            stop 1
         end if
         i = i + 1
         call get_command_argument(i, arg)
         read(arg, *) a
       case default
         print '(a,a,/)', 'Unrecognized command-line option: ', arg
         call print_help()
         stop 1
      end select
      i = i + 1
   end do

   call read_single_column(input, series)

   allocate(result(size(series)))

   call frac_diff_simple(a, series, result)

   call write_table_file(output, result)

   deallocate(series, result)
contains
   subroutine print_help()
      print '(a)', 'Fractionaly differentiate/integrate timeseries'
      print '(a)', ''
      print '(a)', 'usage: fracdiff [options]'
      print '(a)', ''
      print '(a)', 'cmdline options:'
      print '(a)', ''
      print '(a)', '  -v, --version     print version information and exit'
      print '(a)', '  -h, --help        print usage information and exit'
      print '(a)', '  -i, --input       select input file'
      print '(a)', '  -o, --output      select output file'
      print '(a)', '  -a, --order       order of differentiation (negative for integration)'
   end subroutine print_help
end program main_fracdiff
