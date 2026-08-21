program main_psd
   use io, only: read_single_column, write_table_file
   use spectra, only: berg_psd, berg_psd_size
   use version, only: max_ver_len, ver
   use precision, only: wp
   implicit none

   character(len=128) :: arg
   character(len=128) :: input, output
   character(len=max_ver_len) :: v
   real(wp), allocatable :: series(:), f(:), P(:)
   integer :: i, m

   input = 'input.dat'
   output = 'output.dat'
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
       case ('-o', '--output')
          if (i >= command_argument_count()) then
             print '(a)', 'Error: missing value for -o'
            stop 1
         end if
         i = i + 1
         call get_command_argument(i, arg)
         read(arg, *) output
       case ('-m', '--order')
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

   allocate(f(berg_psd_size(size(series))), P(berg_psd_size(size(series))))

   call berg_psd(f, P, series, 1.0_wp, m)

   call write_table_file(output, f, P)

   deallocate(series, f, P)
contains
   subroutine print_help()
      print '(a)', 'Calculate PSD using Burg maximum entropy method'
      print '(a)', ''
      print '(a)', 'usage: psd [options]'
      print '(a)', ''
      print '(a)', 'cmdline options:'
      print '(a)', ''
      print '(a)', '  -v, --version     print version information and exit'
      print '(a)', '  -h, --help        print usage information and exit'
      print '(a)', '  -i, --input       select input file'
       print '(a)', '  -o, --output      select output file'
       print '(a)', '  -m, --order       select berg AR filter length'
   end subroutine print_help
end program main_psd
