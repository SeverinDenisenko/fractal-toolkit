module checks
   use, intrinsic :: iso_fortran_env, only: error_unit
   use stdlib_error, only: stdlib_check => check
   implicit none

   private
   public :: check, last_message

   integer, parameter :: default_err_code = -1
   character(256), save :: last_message = ''

contains
   logical function check(condition, msg, code, ierr) result(failed)
      logical, intent(in) :: condition
      character(*), intent(in), optional :: msg
      integer, intent(in), optional :: code
      integer, intent(out), optional :: ierr

      failed = .false.
      if (.not. present(ierr)) then
         call stdlib_check(condition, msg=msg, code=code)
         return
      end if

      ierr = 0
      if (.not. condition) then
         failed = .true.
         if (present(code)) then
            ierr = code
         else
            ierr = default_err_code
         end if
         if (present(msg)) then
            last_message = msg
            write(error_unit, '(A)') 'check: ' // trim(last_message)
         else
            last_message = ''
         end if
      end if
   end function check
end module checks
