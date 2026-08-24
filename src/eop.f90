module eop
   use precision, only: wp
   use checks, only: check
   implicit none

   private
   public :: eop_label, eop_table
   public :: parse_eop_label, eop_unit_scale
   public :: read_eop, write_eop
   public :: eop_nrow, eop_ncol
   public :: find_eop_column, get_eop_column
   public :: interp_eop_column
   public :: cal_to_mjd, mjd_to_cal

   type :: eop_label
      character(:), allocatable :: name
      character(:), allocatable :: param
      character(:), allocatable :: ref
      character(:), allocatable :: cor
      integer :: num = 0
      integer :: power = 0
      logical :: is_er = .false.
   end type eop_label

   type :: eop_table
      character(:), allocatable :: filename
      type(eop_label), allocatable :: labels(:)
      real(wp), allocatable :: values(:, :)
   end type eop_table

   integer, parameter :: max_record_len = 4096
contains
   pure function eop_unit_scale(lab) result(scale)
      type(eop_label), intent(in) :: lab
      real(wp) :: scale

      scale = 10.0_wp**lab%power
   end function eop_unit_scale

   pure integer function eop_ncol(tab)
      type(eop_table), intent(in) :: tab

      eop_ncol = 0
      if (allocated(tab%values)) eop_ncol = size(tab%values, 1)
   end function eop_ncol

   pure integer function eop_nrow(tab)
      type(eop_table), intent(in) :: tab

      eop_nrow = 0
      if (allocated(tab%values)) eop_nrow = size(tab%values, 2)
   end function eop_nrow

   subroutine parse_eop_label(raw, lab, ierr)
      character(*), intent(in) :: raw
      type(eop_label), intent(out) :: lab
      integer, intent(out), optional :: ierr

      character(len(raw)) :: s
      character(:), allocatable :: toks(:)
      integer :: p, d, ios, nt

      if (present(ierr)) ierr = 0
      s = adjustl(raw)
      lab%name = trim(s)
      lab%power = 0
      lab%num = 0
      lab%is_er = .false.

      p = index(s, '*')
      if (p > 0) then
         read(s(p + 1:), *, iostat=ios) lab%power
         if (check(ios == 0, msg="parse_eop_label: bad unit power in '"//trim(lab%name)//"'", ierr=ierr)) return
         s = s(1:p - 1)
         s = adjustl(s)
      end if

      d = index(s, '.', back=.true.)
      if (d > 1 .and. len_trim(s(d + 1:)) > 0) then
         if (verify(s(d + 1:), '0123456789') == 0) then
            read(s(d + 1:), *) lab%num
            s = s(1:d - 1)
            s = trim(s)
         end if
      end if

      call split_char(trim(s), '_', toks)
      nt = size(toks)

      lab%param = toks(1)
      select case (toks(1))
      case ('COR', 'DA')
         lab%param = trim(s)
      case default
         if (nt >= 2) then
            if (toks(nt) == 'ER' .and. nt >= 2) then
               lab%is_er = .true.
               if (nt >= 3) lab%ref = toks(2)
               if (nt >= 4) lab%cor = toks(3)
            else
               lab%ref = toks(2)
               if (nt >= 3) lab%cor = toks(3)
            end if
         end if
      end select
   end subroutine parse_eop_label

   subroutine find_eop_column(tab, key, idx, ierr)
      type(eop_table), intent(in) :: tab
      character(*), intent(in) :: key
      integer, intent(out) :: idx
      integer, intent(out), optional :: ierr

      integer :: i
      character(len(key)) :: k
      character(:), allocatable :: p

      k = upper(key)
      idx = 0
      if (present(ierr)) ierr = 0
      do i = 1, size(tab%labels)
         p = upper(tab%labels(i)%param)
         if (upper(tab%labels(i)%name) == k .or. p == k) then
            idx = i
            return
         end if
         da_check: if (len(p) > 3) then
            if (p(1:3) == 'DA_' .and. p(4:) == k) then
               idx = i
               return
            end if
         end if da_check
      end do
      if (check(.false., msg="find_eop_column: column '"//trim(key)//"' not found", ierr=ierr)) return
   end subroutine find_eop_column

   subroutine get_eop_column(tab, key, vals, ierr)
      type(eop_table), intent(in) :: tab
      character(*), intent(in) :: key
      real(wp), allocatable, intent(out) :: vals(:)
      integer, intent(out), optional :: ierr

      integer :: idx

      call find_eop_column(tab, key, idx, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      allocate(vals(eop_nrow(tab)))
      vals = tab%values(idx, :)*eop_unit_scale(tab%labels(idx))
   end subroutine get_eop_column

   subroutine interp_eop_column(tab, ix, iy, xq, yq, ierr)
      type(eop_table), intent(in) :: tab
      integer, intent(in) :: ix
      integer, intent(in) :: iy
      real(wp), intent(in) :: xq(:)
      real(wp), intent(out) :: yq(:)
      integer, intent(out), optional :: ierr

      integer :: i, j, n, lo, hi
      real(wp) :: t

      n = eop_nrow(tab)
      if (check(size(xq) == size(yq), msg="interp_eop_column: size mismatch", ierr=ierr)) return
      if (check(ix >= 1 .and. ix <= eop_ncol(tab), msg="interp_eop_column: bad x column index", ierr=ierr)) return
      if (check(iy >= 1 .and. iy <= eop_ncol(tab), msg="interp_eop_column: bad y column index", ierr=ierr)) return
      if (check(n >= 2, msg="interp_eop_column: need at least two rows", ierr=ierr)) return
      do i = 2, n
         if (check(tab%values(ix, i) > tab%values(ix, i - 1), &
                   msg="interp_eop_column: abscissa is not strictly increasing", ierr=ierr)) return
      end do

      do i = 1, size(xq)
         if (check(xq(i) >= tab%values(ix, 1) .and. xq(i) <= tab%values(ix, n), &
                   msg="interp_eop_column: query point out of range", ierr=ierr)) return

         lo = 1
         hi = n
         do while (hi - lo > 1)
            j = (lo + hi)/2
            if (tab%values(ix, j) <= xq(i)) then
               lo = j
            else
               hi = j
            end if
         end do

         t = (xq(i) - tab%values(ix, lo))/(tab%values(ix, hi) - tab%values(ix, lo))
         yq(i) = tab%values(iy, lo)*(1.0_wp - t) + tab%values(iy, hi)*t
      end do
   end subroutine interp_eop_column

   subroutine read_eop(filename, tab, ierr)
      character(*), intent(in) :: filename
      type(eop_table), intent(out) :: tab
      integer, intent(out), optional :: ierr

      character(max_record_len) :: line
      character(:), allocatable :: toks(:)
      integer :: u, ios, i, iline, ncol, nrow
      logical :: is_data

      open(newunit=u, file=filename, status='old', action='read', iostat=ios)
      if (check(ios == 0, msg="read_eop: cannot open '"//filename//"'", ierr=ierr)) return

      ncol = 0
      nrow = 0

      do
         read(u, '(A)', iostat=ios) line
         if (ios /= 0) exit
         if (is_skippable(line)) cycle

         toks = split_ws(line)
         if (size(toks) == 0) cycle

         is_data = token_is_real(toks(1))

         if (.not. is_data .and. ncol == 0) then
            ncol = size(toks)
            allocate(tab%labels(ncol))
            do i = 1, ncol
               call parse_eop_label(toks(i), tab%labels(i), ierr=ierr)
               if (present(ierr) .and. ierr /= 0) then
                  close(u)
                  return
               end if
            end do
         else if (is_data) then
            if (ncol == 0) ncol = size(toks)
            nrow = nrow + 1
         end if
      end do

      if (check(ncol > 0, msg="read_eop: no header or data found in '"//filename//"'", ierr=ierr)) then
         close(u)
         return
      end if
      if (check(nrow > 0, msg="read_eop: no data rows found in '"//filename//"'", ierr=ierr)) then
         close(u)
         return
      end if

      if (.not. allocated(tab%labels)) then
         if (.not. comment_header(u, ncol, tab)) then
            allocate(tab%labels(ncol))
            do i = 1, ncol
               tab%labels(i)%name = 'COL'//int_to_str(i)
               tab%labels(i)%param = tab%labels(i)%name
            end do
         end if
      end if

      allocate(tab%values(ncol, nrow))

      rewind(u)
      iline = 0
      i = 0
      do
         read(u, '(A)', iostat=ios) line
         if (ios /= 0) exit
         iline = iline + 1
         if (is_skippable(line)) cycle

         toks = split_ws(line)
         if (size(toks) == 0) cycle
         if (.not. token_is_real(toks(1))) cycle

         i = i + 1
         read(line, *, iostat=ios) tab%values(:, i)
         if (check(ios == 0, msg="read_eop: malformed data on line "//int_to_str(iline)//" of '"//filename//"'", ierr=ierr)) then
            close(u)
            return
         end if

         if (i == nrow) exit
      end do

      tab%filename = filename
      close(u)
   end subroutine read_eop

   subroutine write_eop(filename, tab, ierr)
      character(*), intent(in) :: filename
      type(eop_table), intent(in) :: tab
      integer, intent(out), optional :: ierr

      integer :: u, ios, i

      open(newunit=u, file=filename, status='replace', action='write', iostat=ios)
      if (check(ios == 0, msg="write_eop: cannot open '"//filename//"'", ierr=ierr)) return

      do i = 1, size(tab%labels)
         if (i < size(tab%labels)) then
            write(u, '(A)', advance='no', iostat=ios) trim(tab%labels(i)%name)//' '
         else
            write(u, '(A)', advance='yes', iostat=ios) trim(tab%labels(i)%name)
         end if
         if (check(ios == 0, msg="write_eop: error writing header", ierr=ierr)) then
            close(u)
            return
         end if
      end do

      do i = 1, eop_nrow(tab)
         write(u, '(*(ES24.16))', iostat=ios) tab%values(:, i)
         if (check(ios == 0, msg="write_eop: error writing row "//int_to_str(i), ierr=ierr)) then
            close(u)
            return
         end if
      end do

      close(u)
   end subroutine write_eop

   pure real(wp) function cal_to_mjd(yr, mo, dy, hr, mn, sec) result(mjd)
      integer, intent(in) :: yr
      integer, intent(in) :: mo
      integer, intent(in) :: dy
      integer, intent(in) :: hr
      integer, intent(in) :: mn
      real(wp), intent(in) :: sec

      integer :: y, m, a, b
      real(wp) :: jd0, frac

      y = yr
      m = mo
      if (m <= 2) then
         y = y - 1
         m = m + 12
      end if

      a = y/100
      b = 2 - a + a/4

      jd0 = int(365.25_wp*real(y + 4716, wp)) + int(30.6001_wp*real(m + 1, wp)) &
            + real(dy, wp) + real(b, wp) - 1524.5_wp

      frac = (real(hr, wp) + real(mn, wp)/60.0_wp + sec/3600.0_wp)/24.0_wp
      mjd = jd0 + frac - 2400000.5_wp
   end function cal_to_mjd

   pure subroutine mjd_to_cal(mjd, yr, mo, dy, hr, mn, sec)
      real(wp), intent(in) :: mjd
      integer, intent(out) :: yr
      integer, intent(out) :: mo
      integer, intent(out) :: dy
      integer, intent(out) :: hr
      integer, intent(out) :: mn
      real(wp), intent(out) :: sec

      real(wp) :: z, f, jd, c, d, e, b, day
      integer :: alpha, a, zint

      jd = mjd + 2400000.5_wp
      z = floor(jd)
      f = jd - z
      zint = int(z)

      if (f >= 0.5_wp) then
         zint = zint + 1
         f = f - 0.5_wp
      else
         f = f + 0.5_wp
      end if

      if (zint < 2299161) then
         a = zint
      else
         alpha = int((real(zint, wp) - 1867216.25_wp)/36524.25_wp)
         a = zint + 1 + alpha - alpha/4
      end if

      b = real(a, wp) + 1524.0_wp
      c = floor((b - 122.1_wp)/365.25_wp)
      d = floor(365.25_wp*c)
      e = floor((b - d)/30.6001_wp)

      day = b - d - floor(30.6001_wp*e) + f
      dy = int(day)
      sec = (day - real(dy, wp))*86400.0_wp

      if (e < 14) then
         mo = int(e - 1.0_wp)
      else
         mo = int(e - 13.0_wp)
      end if

      if (mo > 2) then
         yr = int(c) - 4716
      else
         yr = int(c) - 4715
      end if

      hr = int(sec/3600.0_wp)
      sec = sec - 3600.0_wp*real(hr, wp)
      mn = int(sec/60.0_wp)
      sec = sec - 60.0_wp*real(mn, wp)
   end subroutine mjd_to_cal

   logical function comment_header(u, ncol, tab)
      integer, intent(in) :: u
      integer, intent(in) :: ncol
      type(eop_table), intent(inout) :: tab

      comment_header = comment_header_scan(u, ncol, tab, .true.)
      if (.not. comment_header) comment_header = comment_header_scan(u, ncol, tab, .false.)
   end function comment_header

   logical function comment_header_scan(u, ncol, tab, require_unique)
      integer, intent(in) :: u
      integer, intent(in) :: ncol
      type(eop_table), intent(inout) :: tab
      logical, intent(in) :: require_unique

      character(max_record_len) :: line
      character(:), allocatable :: toks(:)
      integer :: ios, i, j, k, e

      comment_header_scan = .false.
      rewind(u)
      do
         read(u, '(A)', iostat=ios) line
         if (ios /= 0) return
         k = verify(line, ' ')
         if (k == 0) cycle
         if (line(k:k) /= '#') cycle

         toks = split_ws(line(k + 1:))
         if (size(toks) /= ncol) cycle

         do i = 1, ncol
            if (.not. name_like(toks(i))) exit
         end do
         if (i <= ncol) cycle

         if (require_unique) then
            do i = 1, ncol
               do j = i + 1, ncol
                  if (upper(toks(i)) == upper(toks(j))) exit
               end do
               if (j <= ncol) exit
            end do
            if (i <= ncol) cycle
         end if

         allocate(tab%labels(ncol))
         do i = 1, ncol
            e = 0
            call parse_eop_label(toks(i), tab%labels(i), ierr=e)
            if (e /= 0) then
               deallocate(tab%labels)
               return
            end if
         end do
         comment_header_scan = .true.
         return
      end do
   end function comment_header_scan

   pure logical function name_like(tok)
      character(*), intent(in) :: tok

      integer :: i, n
      character :: c
      logical :: has_letter

      name_like = .false.
      n = len_trim(tok)
      if (n == 0 .or. n > 64) return
      if (token_is_real(tok)) return

      has_letter = .false.
      do i = 1, n
         c = tok(i:i)
         select case (c)
         case ('a':'z', 'A':'Z')
            has_letter = .true.
         case ('0':'9', '_', '-', '.', '*')
         case default
            return
         end select
      end do
      name_like = has_letter
   end function name_like

   logical function is_skippable(line)
      character(*), intent(in) :: line

      integer :: k

      k = verify(line, ' ')
      is_skippable = k == 0
      if (.not. is_skippable) is_skippable = line(k:k) == '#'
   end function is_skippable

   pure logical function token_is_real(tok)
      character(*), intent(in) :: tok

      real(wp) :: x
      integer :: ios

      read(tok, *, iostat=ios) x
      token_is_real = ios == 0
   end function token_is_real

   function split_ws(line) result(toks)
      character(*), intent(in) :: line
      character(:), allocatable :: toks(:)

      integer :: i, n, nt, maxl, seglen
      logical :: in_tok

      nt = 0
      maxl = 0
      in_tok = .false.
      do i = 1, len_trim(line)
         if (line(i:i) == ' ' .or. line(i:i) == achar(9)) then
            in_tok = .false.
         else
            if (.not. in_tok) nt = nt + 1
            in_tok = .true.
         end if
      end do

      if (nt == 0) then
         allocate(character(1) :: toks(0))
         return
      end if

      maxl = 1
      n = 0
      in_tok = .false.
      do i = 1, len_trim(line)
         if (line(i:i) == ' ' .or. line(i:i) == achar(9)) then
            if (in_tok) maxl = max(maxl, seglen)
            in_tok = .false.
         else
            if (.not. in_tok) then
               n = n + 1
               seglen = 0
            end if
            seglen = seglen + 1
            in_tok = .true.
         end if
      end do
      if (in_tok) maxl = max(maxl, seglen)

      allocate(character(maxl) :: toks(nt))
      n = 0
      in_tok = .false.
      seglen = 0
      do i = 1, len_trim(line)
         if (line(i:i) == ' ' .or. line(i:i) == achar(9)) then
            in_tok = .false.
         else
            if (.not. in_tok) then
               n = n + 1
               seglen = 0
               toks(n) = ''
            end if
            seglen = seglen + 1
            toks(n)(seglen:seglen) = line(i:i)
            in_tok = .true.
         end if
      end do
   end function split_ws

   subroutine split_char(str, sep, toks)
      character(*), intent(in) :: str
      character, intent(in) :: sep
      character(:), allocatable, intent(out) :: toks(:)

      integer :: i, n, nt, start

      nt = 1
      do i = 1, len_trim(str)
         if (str(i:i) == sep) nt = nt + 1
      end do

      if (len_trim(str) == 0) then
         allocate(character(1) :: toks(1))
         toks(1) = ''
         return
      end if

      allocate(character(len(str)) :: toks(nt))
      toks = ''
      n = 1
      start = 1
      do i = 1, len_trim(str)
         if (str(i:i) == sep) then
            toks(n) = str(start:i - 1)
            n = n + 1
            start = i + 1
         end if
      end do
      toks(n) = str(start:len_trim(str))
      do i = 1, nt
         toks(i) = trim(toks(i))
      end do
   end subroutine split_char

   pure function upper(s) result(u)
      character(*), intent(in) :: s
      character(len(s)) :: u

      integer :: i, ic

      do i = 1, len(s)
         ic = iachar(s(i:i))
         if (ic >= iachar('a') .and. ic <= iachar('z')) then
            u(i:i) = achar(ic - 32)
         else
            u(i:i) = s(i:i)
         end if
      end do
   end function upper

   pure function int_to_str(i) result(s)
      integer, intent(in) :: i
      character(32) :: s

      write(s, '(I0)') i
      s = trim(s)
   end function int_to_str
end module eop
