module io
   use precision, only: wp
   use stdlib_error, only: check
   implicit none

   integer :: default_unit = 42
contains
   subroutine write_table_file(name, array1, array2)
      character(len=*), intent(in) :: name
      real(wp), intent(in) :: array1(:)
      real(wp), optional, intent(in) :: array2(:)

      integer :: i

      if(present(array2)) then
         call check(size(array1) == size(array2), msg="write_table_file: size mismatch")
      end if

      open(unit=default_unit, file=name, status='replace')
      do i = 1, size(array1)
         if(present(array2)) then
            write(default_unit, '(ES20.12, 2X, ES20.12)') array1(i), array2(i)
         else
            write(default_unit, '(ES20.12)') array1(i)
         end if
      end do
      close(default_unit)
   end subroutine write_table_file
end module io
