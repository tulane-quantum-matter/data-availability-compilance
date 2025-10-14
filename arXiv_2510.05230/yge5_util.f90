subroutine tic
   call cpu_time(ticker)
   ticker_set = .true.
end subroutine tic

subroutine toc(message)
   real :: tocker
   character(len=*), optional, intent(in) :: message

   if (ticker_set) then
      call cpu_time(tocker)
      if (present(message)) then
         write (error_unit, *) 'Time lapsed:', tocker - ticker, message
      else
         write (error_unit, *) 'Time lapsed:', tocker - ticker
      end if
   else
      call tic
      write (error_unit, *) 'No ticker set. Now set.'
   end if
end subroutine toc

pure function full_to_sparse(A) result(b)
   implicit none
   complex(8), intent(in) :: A(:, :)
   complex(8), allocatable :: v(:)
   integer :: i
   integer, allocatable :: ii(:)
   type(Spmatz) :: b

   v = reshape(a, [size(a)]); 
   ii = pack([(i, i=1, size(v))], abs(v) > 0); 
   b%nnz = size(ii); 
   allocate (b%entry(b%nnz))
   b%m = size(a, 1)
   b%n = size(a, 2)

   do i = 1, size(ii)
      b%entry(i)%val = v(ii(i))
      b%entry(i)%i = mod(ii(i) - 1, b%m) + 1
      b%entry(i)%j = floor((ii(i) - 1)*1.d0/(b%m)) + 1
   end do
end function

pure function sparse_to_full(b) result(A)
   implicit none
   complex(8), allocatable :: A(:, :)
   integer :: i
   type(Spmatz), intent(in) :: b

   allocate(A(b%m,b%n))
   A = 0
   do i=1,b%nnz
      A(b%entry(i)%i,b%entry(i)%j) = b%entry(i)%val
   end do
end function

subroutine write_spmat_file(spmat, fname)
   use iso_fortran_env
   implicit none
   integer :: unit
   type(Spmatz), intent(in) :: spmat
   character(len=*), intent(in) :: fname

   open (newunit=unit, file=fname, status="replace", position="rewind", action='write')
   call write_spmat_unit(spmat, unit)
   close (unit)
end subroutine

subroutine write_spmat_unit(spmat, unit)
   ! write sparse matrix, compatible with Matrix Market format
   ! import with `Import[filename,"MTX"]` in mathematica
   use iso_fortran_env
   implicit none
   integer, optional, intent(in) :: unit
   type(Spmatz), intent(in) :: spmat
   integer :: i, fu

   if (.not. present(unit)) then
      fu = output_unit
   else
      fu = unit
   end if
   write (fu, '(g0)') "%%MatrixMarket matrix coordinate complex general"
   write (fu, '(3(g0," "))') spmat%m, spmat%n, spmat%nnz
   do i = 1, spmat%nnz
      write (fu, '(2(g0," "),2ES25.16E3)') spmat%entry(i)%i, spmat%entry(i)%j, real(spmat%entry(i)%val), aimag(spmat%entry(i)%val)
   end do
end subroutine

subroutine read_spmat_file(spmat, fname)
   use iso_fortran_env
   implicit none
   integer :: unit
   type(Spmatz), intent(out) :: spmat
   character(len=*), intent(in) :: fname

   open (newunit=unit, file=fname, status="old", action='read')
   call read_spmat_unit(spmat, unit)
   close (unit)
end subroutine

subroutine read_spmat_unit(spmat, unit)
   use iso_fortran_env
   implicit none
   integer, intent(in) :: unit
   type(Spmatz), intent(inout) :: spmat
   character(len=5000) :: typeinfo, sizeinfo
   integer :: i, j, fu, m, n, k
   logical :: issparse
   integer :: matsym ! 0 = general, 1=hermitian, 2=symmetric, 3=skew-symmetric
   real(8) :: re, im

   read (unit, '(a)') typeinfo
   typeinfo = to_lower(typeinfo)
   issparse = index(typeinfo, 'coordinate') > 0
   if (index(typeinfo, ' hermitian') > 0) then
      matsym = 1
   else if (index(typeinfo, ' skew-symmetric') > 0) then
      matsym = 3
   else if (index(typeinfo, ' symmetric') > 0) then
      matsym = 2
   else
      matsym = 0
   end if

   do
      read (unit, '(a)') sizeinfo
      if (sizeinfo(1:1) /= '%') exit
   end do
   if (issparse) then
      read (sizeinfo, *) spmat%m, spmat%n, spmat%nnz
      allocate (spmat%entry(spmat%nnz))
      do i = 1, spmat%nnz
         if (index(typeinfo, 'complex') > 0) then
            read (unit, *) spmat%entry(i)%i, spmat%entry(i)%j, re, im
            spmat%entry(i)%val = cmplx(re, im, 8)
         else
            read (unit, *) spmat%entry(i)%i, spmat%entry(i)%j, re
            spmat%entry(i)%val = re
         end if
      end do
   else
      read (sizeinfo, *) spmat%m, spmat%n
      spmat%nnz = spmat%m*spmat%n
      allocate (spmat%entry(spmat%nnz))
      k = 1
      do i = 1, m
         do j = 1, n
            if (index(typeinfo, 'complex') > 0) then
               read (unit, *) spmat%entry(k)%i, spmat%entry(k)%j, re, im
               spmat%entry(i)%val = cmplx(re, im, 8)
            else
               read (unit, *) spmat%entry(k)%i, spmat%entry(k)%j, re
               spmat%entry(i)%val = re
            end if
            k = k + 1
         end do
      end do
   end if

   spmat = make_general_spmat(spmat, matsym)

end subroutine

pure function make_general_spmat(spmat, matsym) result(gm)
   type(Spmatz) :: gm
   integer, intent(in) :: matsym
   type(Spmatz), intent(in) :: spmat
   type(Entryz) :: ee
   integer :: i, j, m, n, ndiag, nnz

   ndiag = 0
   nnz = spmat%nnz
   do i = 1, spmat%nnz
      if (spmat%entry(i)%i == spmat%entry(i)%j) ndiag = ndiag + 1
   end do

   select case (matsym)
   case (0) ! general, noop
      gm = spmat
   case (1) ! hermitian
      gm%m = spmat%m; gm%n = spmat%n
      gm%nnz = 2*nnz - ndiag
      allocate (gm%entry(gm%nnz))
      j = 0
      do i = 1, nnz
         ee = spmat%entry(i)
         j = j + 1
         gm%entry(j) = ee
         if (ee%i /= ee%j) then
            j = j + 1
            gm%entry(j)%i = ee%j
            gm%entry(j)%j = ee%i
            gm%entry(j)%val = conjg(ee%val)
         end if
      end do
   case (2) ! symmetric
      gm%m = spmat%m; gm%n = spmat%n
      gm%nnz = 2*nnz - ndiag
      allocate (gm%entry(gm%nnz))
      j = 0
      do i = 1, nnz
         ee = spmat%entry(i)
         j = j + 1
         gm%entry(j) = ee
         if (ee%i /= ee%j) then
            j = j + 1
            gm%entry(j)%i = ee%j
            gm%entry(j)%j = ee%i
            gm%entry(j)%val = ee%val
         end if
      end do
   case (3) ! skew-symmetric
      gm%m = spmat%m; gm%n = spmat%n
      gm%nnz = 2*nnz - ndiag
      allocate (gm%entry(gm%nnz))
      j = 0
      do i = 1, nnz
         ee = spmat%entry(i)
         j = j + 1
         gm%entry(j) = ee
         if (ee%i /= ee%j) then
            j = j + 1
            gm%entry(j)%i = ee%j
            gm%entry(j)%j = ee%i
            gm%entry(j)%val = -ee%val
         end if
      end do
   end select

end function

Pure Function to_lower(str) Result(string)

!   ==============================
!   Changes a string to lower case
!   ==============================

   Implicit None
   Character(*), Intent(In) :: str
   Character(LEN(str)) :: string

   Integer :: ic, i

   Character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
   Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

   string = str
   do i = 1, LEN_TRIM(str)
      ic = INDEX(cap, str(i:i))
      if (ic > 0) string(i:i) = low(ic:ic)
   end do

End Function to_lower

Pure Function to_upper(str) Result(string)

!   ==============================
!   Changes a string to upper case
!   ==============================

   Implicit None
   Character(*), Intent(In) :: str
   Character(LEN(str)) :: string

   Integer :: ic, i

   Character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
   Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

!   Capitalize each letter if it is lowecase
   string = str
   do i = 1, LEN_TRIM(str)
      ic = INDEX(low, str(i:i))
      if (ic > 0) string(i:i) = cap(ic:ic)
   end do

End Function to_upper

pure function is_hermitian(mat)
   implicit none
   logical :: is_hermitian
   complex(8), intent(in) :: mat(:, :)
   integer :: i, j, m, n

   is_hermitian = .true.
   outter: do j = 1, size(mat, 2)
      do i = j, size(mat, 1)
         if (abs(mat(i, j) - conjg(mat(j, i))) > 1d-12) then
            is_hermitian = .false.
            return
         end if
      end do
   end do outter
end function

!! sorting utilities
pure recursive subroutine sort_r(mess)
   integer(4), parameter :: bs = 4
   real(bs), intent(inout) :: mess(:)
   real(bs), allocatable :: left(:), right(:)

   include 'yge5-merge-sort.inc' ! integer(4) :: N, ml, mr, i, j, k

end subroutine sort_r

pure recursive subroutine sort_d(mess)
   integer(4), parameter :: bs = 8
   real(bs), intent(inout) :: mess(:)
   real(bs), allocatable :: left(:), right(:)

   include 'yge5-merge-sort.inc' ! integer(4) :: N, ml, mr, i, j, k

end subroutine sort_d

pure recursive subroutine sort_i(mess)
   integer(4), intent(inout) :: mess(:)
   integer(4), allocatable :: left(:), right(:)

   include 'yge5-merge-sort.inc' ! integer(4) :: N, ml, mr, i, j, k

end subroutine sort_i

pure recursive subroutine sort_i8(mess)
   integer(8), intent(inout) :: mess(:)
   integer(8), allocatable :: left(:), right(:)

   include 'yge5-merge-sort.inc' ! integer(4) :: N, ml, mr, i, j, k

end subroutine sort_i8

pure recursive function sortperm_r(mess) result(perm)
   real(4), intent(in) :: mess(:)
   integer, allocatable :: perm(:)

   include 'yge5-merge-sort-perm.inc'
end function sortperm_r

pure recursive function sortperm_d(mess) result(perm)
   real(8), intent(in) :: mess(:)
   integer, allocatable :: perm(:)

   include 'yge5-merge-sort-perm.inc'
end function sortperm_d

pure recursive function sortperm_i(mess) result(perm)
   integer(4), intent(in) :: mess(:)
   integer, allocatable :: perm(:)

   include 'yge5-merge-sort-perm.inc'
end function sortperm_i

pure recursive function sortperm_i8(mess) result(perm)
   integer(8), intent(in) :: mess(:)
   integer, allocatable :: perm(:)

   include 'yge5-merge-sort-perm.inc'
end function sortperm_i8

