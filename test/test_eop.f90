program eop_test
   use precision, only: wp
   use eop, only: eop_label, eop_table, parse_eop_label, read_eop, write_eop, &
      eop_nrow, eop_ncol, find_eop_column, get_eop_column, interp_eop_column, &
      cal_to_mjd, mjd_to_cal
   use stdlib_error, only: check
   implicit none

   type(eop_table) :: tab
   type(eop_table) :: tab2
   type(eop_label) :: scr
   integer :: u, ios, idx
   integer :: yr, mo, dy, hr, mn
   real(wp) :: sec
   real(wp), allocatable :: vals(:)
   real(wp) :: yq(2), xq(2)

   open(newunit=u, file='test_eop_tmp.dat', status='replace', action='write', iostat=ios)
   call check(ios == 0)
   write(u, '(A)') '# sample universal EOP file'
   write(u, '(A)') 'DA_MJD XP YP UT1_UTC LOD XP_ER*-6'
   write(u, '(A)') ''
   write(u, '(A)') '45700.50 -0.132809  0.092060  0.3949652  0.0016989  1368'
   write(u, '(A)') '45701.50 -0.136163  0.094666  0.3933000  0.0016343  1536'
   write(u, '(A)') '45702.50 -0.139500  0.097000  0.3910000  0.0017000  1600'
   close(u)

   call read_eop('test_eop_tmp.dat', tab)
   call check(eop_nrow(tab) == 3)
   call check(eop_ncol(tab) == 6)
   call check(trim(tab%labels(1)%name) == 'DA_MJD')
   call check(trim(tab%labels(1)%param) == 'DA_MJD')
   call check(.not. tab%labels(1)%is_er)

   call parse_eop_label('X_RT_ER*-6', scr)
   call check(trim(scr%name) == 'X_RT_ER*-6')
   call check(trim(scr%param) == 'X')
   call check(trim(scr%ref) == 'RT')
   call check(scr%num == 0)
   call check(scr%power == -6)
   call check(scr%is_er)

   call parse_eop_label('LOD_R.2010', scr)
   call check(trim(scr%name) == 'LOD_R.2010')
   call check(trim(scr%param) == 'LOD')
   call check(scr%num == 2010)
   call check(.not. scr%is_er)

   call parse_eop_label('COR_XP_YP', scr)
   call check(trim(scr%param) == 'COR_XP_YP')

   call find_eop_column(tab, 'ut1_utc', idx)
   call check(idx == 4)
   call find_eop_column(tab, 'MJD', idx)
   call check(idx == 1)

   call get_eop_column(tab, 'XP', vals)
   call check(abs(vals(1) + 0.132809_wp) < 1e-12_wp)
   call check(abs(vals(3) + 0.139500_wp) < 1e-12_wp)

   call get_eop_column(tab, 'XP_ER*-6', vals)
   call check(abs(vals(1) - 1368.0_wp*1e-6_wp) < 1e-12_wp)
   call check(abs(vals(3) - 1600.0_wp*1e-6_wp) < 1e-12_wp)

   xq = [45700.75_wp, 45702.0_wp]
   call interp_eop_column(tab, 1, 2, xq, yq)
   call check(abs(yq(1) - (-0.132809_wp - 0.25_wp*(0.136163_wp - 0.132809_wp))) < 1e-9_wp)
   call check(abs(yq(2) - (-0.136163_wp - 0.5_wp*(0.139500_wp - 0.136163_wp))) < 1e-9_wp)

   call cal_to_mjd_check()

   call read_eop_noheader()
   call read_eop_comment_header()

   call write_eop('test_eop_tmp2.dat', tab)
   call read_eop('test_eop_tmp2.dat', tab2)
   call check(eop_nrow(tab2) == eop_nrow(tab))
   call check(eop_ncol(tab2) == eop_ncol(tab))
   call check(maxval(abs(tab2%values - tab%values)) < 1e-12_wp)

   open(newunit=u, file='test_eop_tmp.dat', status='old')
   close(u, status='delete')
   open(newunit=u, file='test_eop_tmp2.dat', status='old')
   close(u, status='delete')

contains
   subroutine cal_to_mjd_check()
      call check(abs(cal_to_mjd(1984, 1, 1, 12, 0, 0.0_wp) - 45700.50_wp) < 1e-8_wp)
      call check(abs(cal_to_mjd(2000, 1, 1, 12, 0, 0.0_wp) - 51544.50_wp) < 1e-8_wp)
      call check(abs(cal_to_mjd(2026, 8, 24, 0, 0, 0.0_wp) - 61276.0_wp) < 1e-8_wp)

      call mjd_to_cal(45700.50_wp, yr, mo, dy, hr, mn, sec)
      call check(yr == 1984 .and. mo == 1 .and. dy == 1)
      call check(hr == 12 .and. mn == 0 .and. abs(sec) < 1e-6_wp)

      call mjd_to_cal(51544.50_wp, yr, mo, dy, hr, mn, sec)
      call check(yr == 2000 .and. mo == 1 .and. dy == 1)
      call check(hr == 12 .and. mn == 0 .and. abs(sec) < 1e-6_wp)
   end subroutine cal_to_mjd_check

   subroutine read_eop_noheader()
      type(eop_table) :: t3
      real(wp), allocatable :: c(:)

      open(newunit=u, file='test_eop_tmp3.dat', status='replace', action='write')
      write(u, '(A)') '# legacy style file without header line'
      write(u, '(A)') '49187.500 -0.072838  0.233455'
      write(u, '(A)') '49188.500 -0.073601  0.234985'
      close(u)

      call read_eop('test_eop_tmp3.dat', t3)
      call check(eop_nrow(t3) == 2)
      call check(eop_ncol(t3) == 3)
      call check(trim(t3%labels(1)%name) == 'COL1')

      call get_eop_column(t3, 'COL3', c)
      call check(abs(c(2) - 0.234985_wp) < 1e-12_wp)

      open(newunit=u, file='test_eop_tmp3.dat', status='old')
      close(u, status='delete')
   end subroutine read_eop_noheader

   subroutine read_eop_comment_header()
      type(eop_table) :: t4
      real(wp), allocatable :: c(:)

      open(newunit=u, file='test_eop_tmp4.dat', status='replace', action='write')
      write(u, '(A)') '# GPS EARTH ROTATION DATA IN THE IERS FORMAT : EOP(GRGS)'
      write(u, '(A)') '#'
      write(u, '(A)') '# SECONDS SECONDS SECONDS SECONDS OF ARC OF ARC'
      write(u, '(A)') '#  MJD         PM-X      PM-Y       UT1-UTC       DPSI         DEPS'
      write(u, '(A)') '#   (days)      (")       (")        (s)           (")          (")'
      write(u, '(A)') '   56900.500  0.209557  0.337600  -0.3268510      0.000000     0.000000'
      write(u, '(A)') '   56901.500  0.209659  0.336118  -0.3273840      0.000000     0.000000'
      close(u)

      call read_eop('test_eop_tmp4.dat', t4)
      call check(eop_nrow(t4) == 2)
      call check(eop_ncol(t4) == 6)
      call check(trim(t4%labels(1)%name) == 'MJD')
      call check(trim(t4%labels(2)%name) == 'PM-X')
      call check(trim(t4%labels(6)%name) == 'DEPS')

      call get_eop_column(t4, 'PM-Y', c)
      call check(abs(c(1) - 0.337600_wp) < 1e-12_wp)

      open(newunit=u, file='test_eop_tmp4.dat', status='old')
      close(u, status='delete')
   end subroutine read_eop_comment_header
end program eop_test
