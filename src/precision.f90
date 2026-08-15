module precision
   use iso_fortran_env, only: real64, real32
   implicit none

   integer, parameter :: sp = real32
   integer, parameter :: dp = real64
   integer, parameter :: wp = dp

end module precision
