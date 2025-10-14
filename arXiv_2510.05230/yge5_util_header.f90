   !!!!!!!!!!!!!!!!!!!!
   ! yge5 utilities

   ! timer
   real :: ticker
   logical :: ticker_set = .false.

   ! sparse matrix
   type Entryz
      integer :: i, j
      complex(8) :: val
   end type Entryz

   type Spmatz
      integer :: nnz ! # of nonzero
      integer :: m, n
      type(Entryz), allocatable :: entry(:)
   end type Spmatz

   ! in place sort procedures
   ! using merge sort
   interface sort
      module procedure sort_r, sort_d, sort_i, sort_i8
   end interface sort

   interface sortperm
      ! returns integer (:)
      ! mach dependent kind
      module procedure sortperm_r, sortperm_d, sortperm_i, sortperm_i8
   end interface sortperm

   interface write_spmat
      module procedure write_spmat_unit, write_spmat_file
   end interface write_spmat

   interface read_spmat
      module procedure read_spmat_unit, read_spmat_file
   end interface read_spmat
   !!!!!!!!!!!!!!!!!!!!