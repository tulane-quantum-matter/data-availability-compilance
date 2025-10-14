submodule(Hamiltonian_main) ham_TSC_smod

   Use Operator_mod
   Use WaveFunction_mod
   Use Lattices_v3
   Use MyMats
   Use Random_Wrap
   Use Files_mod
   Use Matrix
   !use control
   Use Observables
   Use Fields_mod
   !Use Predefined_Hoppings
   Use LRC_Mod
   use runtime_error_mod
   use iso_fortran_env

   Implicit none

   integer, parameter :: rk = kind(0.d0), ck = kind(0.d0)
   complex(ck), parameter :: im = (0.0, 1.0) ! imaginary unit

   type, extends(ham_base) :: ham_TSC
   contains
      ! Set Hamiltonian-specific procedures
      procedure, nopass :: Ham_Set
      procedure, nopass :: Alloc_obs
      procedure, nopass :: Obser
      procedure, nopass :: ObserT
      procedure, nopass :: weight_reconstruction
      procedure, nopass :: GR_reconstruction
      procedure, nopass :: GRT_reconstruction
#ifdef HDF5
      procedure, nopass :: write_parameters_hdf5
#endif
   end type ham_TSC

   !#PARAMETERS START# VAR_lattice
   Character(len=64) :: Model = ''  ! Value irrelevant
   Character(len=64) :: Lattice_type = 'Bilayer_square'  ! Possible Values: 'Bilayer_honeycomb' 'Bilayer_square'
   Integer :: L1 = 6   ! Length in direction a_1
   Integer :: L2 = 6   ! Length in direction a_2
   !#PARAMETERS END#

   !#PARAMETERS START# VAR_TSC
   real(Kind=Kind(0.d0)) :: ham_T = 1.d0          ! Hopping parameter
   real(Kind=Kind(0.d0)) :: ham_delta = 0.4       ! SOC parameter
   real(Kind=Kind(0.d0)) :: Ham_chem = -0.5       ! Chemical potential
   real(Kind=Kind(0.d0)) :: Ham_U = 1.7           ! Attractive Hubbard interaction strength (>0)
   real(Kind=Kind(0.d0)) :: Ham_dUedge = 0.d0     ! Hubbard interaction on the edge
   real(Kind=Kind(0.d0)) :: Dtau = 0.05           ! Thereby Ltrot=Beta/dtau
   real(Kind=Kind(0.d0)) :: Beta = 60             ! Inverse temperature
   !logical :: Projector = .false.                ! Whether the projective algorithm is used, DO NOT UNCOMMENT THIS LINE
   real(Kind=Kind(0.d0)) :: Theta = 100           ! Projection parameter
   !logical :: Symm = .false.                      ! Whether symmetrization takes place, DO NOT UNCOMMENT THIS LINE
   Integer :: N_part = -1                         ! Number of particles in trial wave function. If N_part < 0 -> N_part = L1*L2/2
   logical :: Hubbard_Mz = .true.                 ! whether to use Mz HS-decoupling scheme, which may be more stable
   logical :: open_boundary = .false.             ! whether open-boundary condition along zig-zag edge, at y=0-
   logical :: no_self_opx = .false.               ! whether hop or interacts with itself across a periodic boundary
   logical :: no_accumulate_opx = .false.         ! whether setting kinetic and interaction terms accumulatively
   integer :: excluded_sites_input(3000) = 0      ! orbitals to exclude from input file, [1,Ndim]
   integer :: cheatcode(64) = 0                   ! debug code array
   !#PARAMETERS END#

   !!!!!!!!!!!!!!!!!!!!
   ! yge5 utilities
   include 'yge5_util_header.f90'
   !!!!!!!!!!!!!!!!!!!!

   Type(Lattice), target :: Latt
   Type(Unit_cell), target :: Latt_unit
   INTEGER :: nf_calc, nf_reconst

   Integer, allocatable, target :: List(:, :), Invlist(:, :)  ! For orbital structure of Unit cell
   integer, pointer :: unit_cell_list(:, :), unit_cell_inv_list(:, :)

   !!!!!!!!!!!!!!
   ! additional lattice book keeping variables ----
   ! open boundary markers
   integer, allocatable :: edge_cell_list(:), margin_cell_list(:), edge_site_list(:), bottom_row_cell_list(:), top_row_cell_list(:) ! unit cell across the edge
   integer :: n_edge_cell = 0, n_margin_cell = 0, n_edge_site = 0, n_bottom_row_cell, n_top_row_cell
   logical, allocatable :: is_edge_cell(:), is_margin_cell(:), is_edge_site(:), is_excluded_site(:), is_bottom_row(:), is_top_row(:)

   ! excluded site markers
   integer :: n_excluded = 0
   integer, allocatable :: excluded_sites(:)

   real(rk), allocatable, target :: xys(:, :) ! lattice coordinates
   real(rk), pointer :: xs(:), ys(:)

   ! 1-d list of orders of the 1..ndim orbitals by (1..Norb, 1..L1, 1..L2)
   integer, allocatable :: ordered_site_list(:)
   ! 1-d list of orders of the particle layer orbitals by (1..Nsublattice, 1..L1, 1..L2)
   integer, allocatable :: ordered_orbital_list(:)

   !!!!!!!!
   ! Hamiltonian terms matrices

   ! noninteracting terms, book keeping
   type(Spmatz), allocatable :: hopmat(:), pairmat(:)

   ! interaction terms, book keeping
   integer :: N_Hubbard = 0, N_Hubbard_edge = 0

   !!!!!!!!
   ! Observables

   ! --- For obs_scal ---
   enum, bind(c)
      enumerator :: os_kin = 1, os_pot, os_part, os_ener, os_hop &
         , os_pair, os_hubd, os_deltaS2, os_gr12, os_deltaS12 &
         , os_deltaS1h, os_deltaS2h, os_deltaS3h, os_deltaS4h &
         , os_spinxh, os_spinyh, os_spinzh, os_denh &
         , os_mstat, os_mstat2, os_mall &
         , os_grall, os_deltaSall
   end enum
   integer :: os_last = os_deltaSall ! The last scalar observable to output

   ! --- For obs_eq ---
   enum, bind(c)
      enumerator :: oe_deltas_intra = 1, oe_deltas_inter, oe_gr_intra, oe_gr_inter
   end enum
   integer :: oe_last = oe_gr_inter ! The last obs_eq to output
   integer, allocatable :: interYs(:)

   ! --- For obs_tau ---
   enum, bind(c)
      enumerator :: ot_gr_intrawire = 1, ot_gr_local_edge, ot_gr_local_bulk
   end enum
   integer :: ot_last = ot_gr_local_bulk ! The last obs_tau to output

contains

   module Subroutine Ham_Alloc_TSC
      allocate (ham_TSC :: ham)
   end Subroutine Ham_Alloc_TSC

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_TSC_read_write_parameters.F90"

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the Hamiltonian
!--------------------------------------------------------------------
   Subroutine Ham_Set

#if defined (MPI) || defined(TEMPERING)
      Use mpi
#endif
      Implicit none

      integer :: ierr, nf, unit_info
      Character(len=64) :: file_info

#ifdef MPI
      Integer :: Isize, Irank, irank_g, isize_g, igroup
      Integer :: STATUS(MPI_STATUS_SIZE)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
      call MPI_Comm_rank(Group_Comm, irank_g, ierr)
      call MPI_Comm_size(Group_Comm, isize_g, ierr)
      igroup = irank/isize_g
#endif

      call tic

      ! From dynamically generated file "Hamiltonian_TSC_read_write_parameters.F90"
      call read_parameters()

      Ltrot = nint(beta/dtau)
      Thtrot = 0
      if (Projector) Thtrot = nint(theta/dtau)
      Ltrot = Ltrot + 2*Thtrot

      ! Setup the Bravais lattice
      Call Ham_Latt

      ! always half-filled in a nambu setting
      if (N_part < 0) N_part = Ndim/2
      N_SUN = 1
      N_FL = 2

      ! Setup the hopping / single-particle part
      Call Ham_Hop

      ! Setup the interaction.
      call Ham_V_builder
      write (error_unit, *) 'Cont after Ham_V_builder'

      ! Setup the trival wave function, in case of a projector approach
      if (Projector) then
         Call Ham_Trial()
         write (error_unit, *) 'Trial wavefunction generated'
      end if
#ifdef MPI
      If (Irank_g == 0) then
#endif
         File_info = "info"
#if defined(TEMPERING)
         write (File_info, '(A,I0,A)') "Temp_", igroup, "/info"
#endif

         Open (newunit=unit_info, file=file_info, status="unknown", position="append")
         Write (unit_info, *) '====================================='
         Write (unit_info, *) 'Model is      : TSC'
         Write (unit_info, *) 'Lattice is    : ', Lattice_type
         Write (unit_info, *) 'L1            : ', L1
         Write (unit_info, *) 'L2            : ', L2
         Write (unit_info, *) '# of orbitals : ', Ndim
         Write (unit_info, *) '# of flavors  : ', N_FL
         Write (unit_info, *) 'Symm. decomp  : ', Symm
         if (Projector) then
            Write (unit_info, *) 'Projective version'
            Write (unit_info, *) 'Theta         : ', Theta
            Write (unit_info, *) 'Tau_max       : ', beta
            Write (unit_info, *) '# of particles: ', N_part
         else
            Write (unit_info, *) 'Finite temperture version'
            Write (unit_info, *) 'Beta          : ', Beta
         end if
         Write (unit_info, *) 'dtau,Ltrot_eff: ', dtau, Ltrot
         Write (unit_info, *) 't             : ', Ham_T
         Write (unit_info, *) 'delta         : ', ham_delta
         Write (unit_info, *) 'Ham_U         : ', Ham_U
         Write (unit_info, *) 'Ham_dUedge    : ', Ham_dUedge
         Write (unit_info, *) 'Ham_chem      : ', Ham_chem
         write (unit_info, *) 'Hubbard_Mz    : ', Hubbard_Mz

         if (.not. Hubbard_Mz) then
    write (error_unit, *) 'Using SU(2)-symmetric Hubbard-decoupling for testing purpose. Loss of Green function precision expected.'
         end if

         if (open_boundary) then
            write (unit_info, *) 'Cylinder geometry with open boundary at y=0.'
         else
            write (unit_info, *) 'Torus geometry, plainly periodic in both directions.'
         end if
         if (no_accumulate_opx) then
            write (unit_info, *) 'Do not accumulate hamiltonian couplings over the same link/orbital sets.'
         else
            write (unit_info, *) 'Allow to accumulate hamiltonian couplings over the same link/orbital sets.'
         end if
         if (no_self_opx) then
            write (unit_info, *) 'Do not allow coupling to oneself, e.g. across the periodic boundary.'
         else
            write (unit_info, *) 'Allow coupling to oneself, e.g. across the periodic boundary.'
         end if

         write (unit_info, *) 'Number of sites excluded : ', n_excluded

         ! symmetric reconstructions
         if (Ham_U < 0.d0) then
            Write (unit_info, *) 'No proper symmetry assumed for U=', ham_U
            Call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
         end if
         Write (unit_info, *) 'Assuming particle hole symmetry. Sign free Problem.'

         if (Projector) then
            Do nf = 1, N_FL
               Write (unit_info, *) 'Degen of right trial wave function: ', WF_R(nf)%Degen
               Write (unit_info, *) 'Degen of left  trial wave function: ', WF_L(nf)%Degen
            end do
         end if
         Close (unit_info)
#ifdef MPI
      End if
#endif

      ! weight reconstruction setup
      allocate (Calc_Fl(N_FL))
      nf_calc = 1
      nf_reconst = 2
      Calc_Fl(nf_calc) = .True.
      Calc_Fl(nf_reconst) = .False.

   End Subroutine Ham_Set

!--------------------------------------------------------------------
!> Sets  the  Lattice
!--------------------------------------------------------------------
   Subroutine Ham_Latt
      Use Predefined_Lattices

      Implicit none
      integer, allocatable :: il(:)
      integer :: i, j, k, l, m, n, k1, k2, p, q, qs(1)

      write (error_unit, *) 'Building lattice'

      If (Lattice_Type /= "Bilayer_square") then
         Write (output_unit, *) 'The model is only defined for the bilayer square lattice'
         Call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
      End if

      call Predefined_Latt(Lattice_type, L1, L2, Ndim, List, Invlist, Latt, Latt_Unit)
      unit_cell_list => list
      unit_cell_inv_list => invlist

      call fill_latt_invlist

      allocate (xys(ndim, 2))
      do i = 1, Ndim
         n = list(i, 1)
         m = list(i, 2)
         xys(i, :) = latt%list(n, 1)*latt%a1_p + latt%list(n, 2)*latt%a2_p + Latt_unit%Orb_pos_p(m, 1:2)
      end do
      xs => xys(:, 1)
      ys => xys(:, 2)

      call fill_order_lists

      call mark_open_boundary ! mark open boundary across y=0

      call mark_excluded_sites

      write (error_unit, *) 'Lattice generated'

      Open (newunit=j, file='lattice_info', status="unknown", action="write")
      do i = 1, Ndim
         k = list(i, 1)
         n = list(i, 2)
         m = 0 ! edge or marginal cell?
         if (is_bottom_row(k)) then
            m = 2
         elseif (is_top_row(k)) then
            m = 1
         end if
         l = 0
         if (is_edge_site(i)) l = 1 ! edge site?
         p = 0
         if (is_excluded_site(i)) p = 1 ! excluded?
         write (j, "(10(g0,' '))") i, k, n &   ! orbital, unit cell, suborbital indices
            , latt%list(k, 1), latt%list(k, 2) &   ! lattice vector coordinates
            , xs(i), ys(i) &   ! x, y
            , m, l, p   ! upper/lower edge cell, edge site, excluded site
      end do
      close (j)

      Open (newunit=j, file='lattice_NN_info', status="unknown", action="write")
      do i = 1, Ndim
         k = list(i, 1)
         n = list(i, 2)

         write (j, "(14(g0,' '))") i, k, n, & ! orbital, unit cell, suborbital indices
            latt%list(k, 1), latt%list(k, 2), & ! lattice vector coordinates
            [(latt%nnlist(k, -1, l), l=-1, 1)], &
            [(latt%nnlist(k, 0, l), l=-1, 1)], &
            [(latt%nnlist(k, 1, l), l=-1, 1)]
      end do
      close (j)

      Open (newunit=j, file='lattice_unit_info', status="unknown", action="write")
      do i = 1, Ndim
         k = list(i, 1)
         n = list(i, 2)

         write (j, "(14(g0,' '))") i, k, n, & ! orbital, unit cell, suborbital indices
            latt%list(k, 1), latt%list(k, 2), & ! lattice vector coordinates
            [(latt%nnlist(k, -1, l), l=-1, 1)], &
            [(latt%nnlist(k, 0, l), l=-1, 1)], &
            [(latt%nnlist(k, 1, l), l=-1, 1)]
      end do
      close (j)

      open (newunit=j, file='lattice_k_info', status='unknown', action='write')
      do i = 1, latt%n
         write (j, *) modulo(latt%listk(i, 1), l1), modulo(latt%listk(i, 2), l2)
      end do
      close (j)

   contains

      ! fill the ordered_*_list in orbital_1, x_1=0,..,L1-1, y_1=0,..L2-1 then orbital_2, x_2, y_2
      subroutine fill_order_lists
         integer :: i, j, k

         ! use radix sort to sort the indexes
         ordered_orbital_list = sortperm(modulo(latt%list(list([(i, i=1, ndim)], 1), 1), L1))
         ordered_orbital_list = ordered_orbital_list( &
                                sortperm(modulo(latt%list(list(ordered_orbital_list, 1), 2), L2)))
         ordered_site_list = pack(ordered_orbital_list, &
                                  modulo(ordered_orbital_list - 1, latt_unit%Norb) <= latt_unit%Norb/2 - 1)

         if (cheatcode(2) /= 0) then
            Open (newunit=i, file='lattice_ordered_orbital_list', status="unknown", action="write")
            do k = 1, size(ordered_orbital_list)
               j = ordered_orbital_list(k)
               write (i, '(4(g0,1X))') j, list(j, 2), latt%list(list(j, 1), 1), latt%list(list(j, 1), 2)
            end do
            close (i)
            Open (newunit=i, file='lattice_ordered_site_list', status="unknown", action="write")
            do k = 1, size(ordered_site_list)
               j = ordered_site_list(k)
               write (i, '(3(g0,1X))') j, latt%list(list(j, 1), 1), latt%list(list(j, 1), 2)
            end do
            close (i)
         end if

      end subroutine fill_order_lists

      subroutine mark_open_boundary
         write (error_unit, *) 'Labeling open boundary'

         allocate (is_edge_cell(latt%N), is_margin_cell(latt%N), is_bottom_row(latt%N), is_top_row(latt%N))
         is_edge_cell = .false.
         is_margin_cell = .false.
         allocate (is_edge_site(Ndim))
         is_edge_site = .false.
         n_edge_cell = 0
         n_margin_cell = 0
         if (open_boundary) then
            il = [(i, i=1, latt%N)]
            is_bottom_row = (modulo(latt%list(:, 2), L2) == 0)
            is_top_row = modulo(latt%list(:, 2), L2) == modulo(-1, L2)
            is_edge_cell = is_bottom_row .or. is_top_row
            edge_cell_list = pack(il, is_edge_cell)
            bottom_row_cell_list = pack(il, is_bottom_row)
            top_row_cell_list = pack(il, is_top_row)
            n_edge_cell = size(edge_cell_list, 1)
            n_bottom_row_cell = size(bottom_row_cell_list, 1)
            n_top_row_cell = size(top_row_cell_list, 1)

            is_margin_cell = is_edge_cell

            do i = 1, n_edge_cell
               j = edge_cell_list(i)
               do k = 1, Latt_unit%Norb
                  is_edge_site(Invlist(j, k)) = .true.
               end do
            end do
            margin_cell_list = pack(il, is_margin_cell)
            n_margin_cell = size(margin_cell_list)
            edge_site_list = pack([(i, i=1, Ndim)], is_edge_site)
            n_edge_site = size(edge_site_list)
         end if
      end subroutine mark_open_boundary

      subroutine mark_excluded_sites
         allocate (is_excluded_site(ndim))
         is_excluded_site = .false.
         do i = 1, size(excluded_sites_input, 1)
            if (excluded_sites_input(i) >= 1 .and. excluded_sites_input(i) <= Ndim) then
               is_excluded_site(excluded_sites_input(i)) = .true.
            else
               exit
            end if
         end do
         j = 0
         do i = 1, Ndim ! also set the nambu partners
            if (is_excluded_site(i)) then
               if (latt_unit%norb/i >= 2) then ! particle sector
                  if (.not. (is_excluded_site(i + latt_unit%norb/2))) then
                     j = 1
                     is_excluded_site(i + latt_unit%norb/2) = .true.
                  end if
               else ! hole sector
                  if (.not. (is_excluded_site(i - latt_unit%norb/2))) then
                     j = 1
                     is_excluded_site(i - latt_unit%norb/2) = .true.
                  end if
               end if
            end if
         end do
         n_excluded = count(is_excluded_site)
         excluded_sites = pack([(i, i=1, Ndim)], is_excluded_site)
        if (j > 0) write (error_unit, *) 'Some excluded sites have their Nambu partners unset! We have set them to be excluded too.'
         write (error_unit, *) 'Number of excluded sites: ', n_excluded
      end subroutine mark_excluded_sites

   end Subroutine Ham_Latt

   subroutine fill_latt_invlist
      ! Fill rest of latt%invlist for easy traversal in lattice
      ! By default, latt%invlist(-L1..L1, -L2..L2) only have latt%N (arbitary) nonzero entries
      implicit none
      integer :: i, j, k, m, n

      do i = 1, latt%n
         j = latt%list(i, 1)
         k = latt%list(i, 2)

         do m = j - 3*L1, j + 3*L1, L1
            do n = k - 3*L2, k + 3*L2, L2
               if (m >= -L1 .and. m <= L1 .and. n >= -L2 .and. n <= L2) then
                  if (m /= j .or. n /= k) latt%invlist(m, n) = i
               end if
            end do
         end do
      end do

   end subroutine fill_latt_invlist

!--------------------------------------------------------------------
!> Sets  the Hopping
!--------------------------------------------------------------------

   Subroutine Ham_Hop
      Implicit none
      Integer :: nf, I, Ix, Iy, lmax, n1, n2, n3, n4, j, k, l, m, n
      integer :: right_uc, upper_uc
      real(kind=kind(0.0d0)) :: X, g ! hopping strengths of the NN hopping and the pairing amplitude
      complex(ck) :: t, delta, mu
      complex(ck), allocatable :: opthop(:, :, :), optpair(:, :, :)

      ! no need to divide by two
      t = -ham_T; 
      delta = ham_delta; 
      mu = ham_chem; 
      allocate (opthop(ndim, ndim, n_fl), optpair(ndim, ndim, n_fl))
      opthop = 0
      optpair = 0

      !!!!! Add kinetic terms. Excluded sites removed at the end
      !!!!! nearest neighbor hopping -----

      do nf = 1, N_FL
         Do i = 1, Latt%N
            right_uc = Latt%nnlist(i, 1, 0) ! right unit cell
            upper_uc = Latt%nnlist(i, 0, 1) ! upper unit cell

            ! to right
            ! hop to the right unit
            if ((.not. no_self_opx) .or. right_uc /= i) then
               call set_nn_hopping(opthop(:, :, nf), Invlist(i, 1), Invlist(right_uc, 1), t, no_accumulate_opx)
            end if

            ! to upper
            ! hop to the upper unit
            if ((.not. no_self_opx) .or. upper_uc /= i) then
               if (.not. is_bottom_row(upper_uc)) then
                  call set_nn_hopping(opthop(:, :, nf), Invlist(upper_uc, 1), Invlist(i, 1), t, no_accumulate_opx)
               end if
            end if

            ! chemical potential
            n1 = Invlist(i, 1); n2 = Invlist(i, 2)
            if (no_accumulate_opx) then
               opthop(n1, n1, nf) = 0
               opthop(n2, n2, nf) = 0
            end if
            opthop(n1, n1, nf) = opthop(n1, n1, nf) - mu
            opthop(n2, n2, nf) = opthop(n2, n2, nf) + mu

            ! pairing right, px
            if ((.not. no_self_opx) .or. right_uc /= i) then
               call set_nn_pairing(optpair(:, :, nf), Invlist(i, 1), Invlist(right_uc, 1), delta, no_accumulate_opx)
            end if

            ! pairing upper, i py
            if ((.not. no_self_opx) .or. upper_uc /= i) then
               if (.not. is_bottom_row(upper_uc)) then
           call set_nn_pairing(optpair(:, :, nf), Invlist(upper_uc, 1), Invlist(i, 1), ((-1)**(nf + 1))*im*delta, no_accumulate_opx)
               end if
            end if

         end do
      end do

      !! remove excluded sites ---
      do i = 1, n_excluded
         k = excluded_sites(i)
         opthop(:, k, :) = 0
         opthop(k, :, :) = 0

         optpair(:, k, :) = 0
         optpair(k, :, :) = 0
      end do

      Allocate (Op_T(1, N_FL))  !there are two hopping matrix, the NN hopping and the p wave pairing
      do nf = 1, N_FL
         Call Op_make(Op_T(1, nf), Ndim)
         Op_T(1, nf)%P = [(i, i=1, Ndim)]
         Op_T(1, nf)%g = -Dtau
         op_t(1, nf)%O = opthop(:, :, nf) + optpair(:, :, nf)
      end do

      ! convert to sparse matrix spmatz to consolidate entries ---
      ! used to compute kinetic energies
      allocate (hopmat(n_fl))
      do i = 1, n_fl
         hopmat(i) = full_to_sparse(opthop(:, :, i))
      end do

      allocate (pairmat(n_fl))
      do i = 1, n_fl
         pairmat(i) = full_to_sparse(optpair(:, :, i))
      end do

      ! Write Op_T%O to plain text file, eg for debugging
      call write_opt_mat

      ! check hermicity
      do j = 1, size(Op_T, 1)
         do nf = 1, n_fl
            if (.not. is_hermitian(Op_T(j, nf)%O)) then
               write (error_unit, *) 'Op_T(', j, ',', nf, ')%O is not Hermitian. Aborting.'
               stop
            end if
         end do
      end do

      !! Set the kinetic operators
      do i = 1, size(Op_T, 1)
         do j = 1, size(Op_T, 2)
            Call Op_set(Op_T(i, j))
         end do
      end do
   contains

      subroutine set_nn_hopping(opto, n1, n2, t, no_accumulate_opx)
         ! Add hopping t from A site n1, to B site n2
         integer, intent(in) :: n1, n2
         complex(ck), intent(in) :: t
         logical, intent(in) :: no_accumulate_opx
         complex(ck), intent(inout) :: opto(:, :)
         integer :: u

         u = Latt_unit%norb/2

         if (no_accumulate_opx) opto(n1, n2) = 0
         opto(n1, n2) = opto(n1, n2) + t

         if (no_accumulate_opx) opto(n2, n1) = 0
         opto(n2, n1) = opto(n2, n1) + conjg(t) ! h.c.

         ! nambus
         if (no_accumulate_opx) opto(n1 + u, n2 + u) = 0
         opto(n1 + u, n2 + u) = opto(n1 + u, n2 + u) - conjg(t)

         if (no_accumulate_opx) opto(n2 + u, n1 + u) = 0
         opto(n2 + u, n1 + u) = opto(n2 + u, n1 + u) - t ! h.c.
      end subroutine set_nn_hopping

      subroutine set_nn_pairing(opto, n1, n2, delta, no_accumulate_opx)
         integer, intent(in) :: n1, n2
         complex(ck), intent(in) :: delta
         logical, intent(in) :: no_accumulate_opx
         complex(ck), intent(inout) :: opto(:, :)
         integer :: u

         u = Latt_unit%norb/2

         if (no_accumulate_opx) then
            opto(n1, n2 + u) = 0
            opto(n1 + u, n2) = 0
            opto(n2 + u, n1) = 0
            opto(n2, n1 + u) = 0
         end if
         opto(n1, n2 + u) = -delta
         opto(n2, n1 + u) = delta

         opto(n1 + u, n2) = conjg(delta)
         opto(n2 + u, n1) = -conjg(delta)
      end subroutine set_nn_pairing

      subroutine write_opt_mat
         integer :: i

         ! Write Op_T%O as stored in memory
         Open (newunit=i, file="hopmat.txt", status="unknown", position="rewind")
         call write_spmat(hopmat(1), i)
         close (i)
         Open (newunit=i, file="pairmat.txt", status="unknown", position="rewind")
         call write_spmat(pairmat(1), i)
         close (i)

         ! Write Op_T%O as ordered in physical coordinates
         ! (orb1, x1, y1, orb2, x2, y2)
         Open (newunit=i, file="hopmat-ordered.txt", status="unknown", position="rewind")
         call write_spmat(full_to_sparse(opthop(ordered_orbital_list, ordered_orbital_list, 1)), i)
         close (i)
         Open (newunit=i, file="pairmat-ordered.txt", status="unknown", position="rewind")
         call write_spmat(full_to_sparse(optpair(ordered_orbital_list, ordered_orbital_list, 1)), i)
         close (i)

      end subroutine write_opt_mat

   end Subroutine Ham_Hop

!--------------------------------------------------------------------
!> Sets the trial wave function
!--------------------------------------------------------------------
   Subroutine Ham_Trial()

      ! we generally opt for finite-T
      !Use Predefined_Trial
      use MyMats

      Implicit none
      integer :: s, x, n, i
      complex(ck), allocatable :: U(:, :), H(:, :)
      real(rk), allocatable :: W(:)
      real(rk) :: r(2) = [real(rk) :: 0.d0, 1.d0], v

      allocate (W(Ndim), H(Ndim, Ndim), U(Ndim, Ndim))
      allocate (wf_l(n_fl), wf_r(n_fl))

      Do s = 1, N_fl
         H = 0
         Do i = 1, size(op_t, 1)
            H = H + op_t(i, s)%O/2 ! divide by two for nambu
         end do

         ! random potential to gap and break translational symmetry
         do i = 1, size(H, 1), 2
            CALL RANDOM_NUMBER(r) ! add some random staggerring to lift degeneracy
            v = max(ham_T, ham_U/4)*sign(0.5 + 0.4*abs(r(2) - r(1)), r(2) - r(1))
            v = v*(-1)**(latt%list(list(i, 1), 1) + latt%list(list(i, 1), 2))
            !v = ham_t*0.4*(-1)**(latt%list(list(i, 1), 1) + latt%list(list(i, 1), 2))
            !write (*, *) v
            H(i, i) = H(i, i) + v
            H(i + 1, i + 1) = H(i + 1, i + 1) - v
         end do

         call diag(H, U, W)

         WF_L(s)%P = (U(:, 1:N_part))
         WF_L(s)%Degen = w(N_part + 1) - w(N_part)
         WF_R(s)%P = (conjg(WF_L(s)%P))
         WF_R(s)%Degen = w(N_part + 1) - w(N_part)
      End do

   contains

      pure function arb_permute(U) result(ww)
         complex(ck), intent(in) :: U(:, :)
         complex(ck), allocatable :: ww(:, :)
         integer, parameter :: nn = 4
         integer :: i
         integer, allocatable :: ss(:)

         ss = mod([(list(i, 1), i=1, size(U, 1))], nn); 
         ww = U(sortperm(ss), :)

      end function arb_permute

   end Subroutine Ham_Trial

!--------------------------------------------------------------------
!> Sets the interaction. Hubbard, Hubbard on edge.
!--------------------------------------------------------------------
   Subroutine Ham_V_builder

      !Use Predefined_Int
      Implicit none

      Integer :: nf, I, i2, j, nb
      Real(ck) :: X, U ! for flipping sign between flavors
      complex(ck) :: g ! Strength
      integer :: N_V_total

      call hubbard_V(counter=N_Hubbard, countOnly=.true.)

      N_V_total = N_hubbard
      write (*, *) 'N_V_total = N_hubbard', N_V_total, N_hubbard
      Allocate (Op_V(N_V_total, N_FL))

      call hubbard_V(Op_V(1:N_Hubbard, :), N_Hubbard, countOnly=.false.)

      do j = 1, size(Op_V, 1)
         do nf = 1, n_fl
            if (.not. is_hermitian(Op_V(j, nf)%O)) then
               write (error_unit, *) 'Op_V(', j, ',', nf, ')%O is not Hermitian. Aborting.'
               stop
            end if
            ! write (*,*) op_v(j,nf)%p
            ! write (*,*) op_v(j,nf)%O
         end do
      end do

      write (error_unit, *) 'All interactions generated.'

   end Subroutine Ham_V_builder

!--------------------------------------------------------------------
!> Hubbard interaction. The input has ham_U = -U.
!--------------------------------------------------------------------
   subroutine hubbard_V(Op_V, counter, countOnly)
      logical, intent(in) :: countOnly
      Type(Operator), dimension(:, :), intent(inout), optional :: Op_V
      integer, intent(out) :: counter
      integer :: nf, i, j, k, l
      complex(ck) :: U, g

      counter = 0
      ! attractive Hubbard interaction in the bulk
      if (ham_U /= 0 .or. ham_dUedge /= 0) then
         Do nf = 1, N_FL
            counter = 0
            Do i = 1, Latt%N
               if (.not. is_excluded_site(invlist(i, 1))) then
                  counter = counter + 1
                  j = counter
                  if (.not. countOnly) then
                     Call Op_make(Op_V(j, nf), 2)
                     Op_V(j, nf)%P(1) = Unit_cell_inv_list(i, 1)
                     Op_V(j, nf)%P(2) = Unit_cell_inv_list(i, 2)
                     Op_V(j, nf)%O(1, 1) = 1
                     Op_V(j, nf)%O(2, 2) = -1
                     U = ham_U
                     if (is_edge_site(invlist(i, 1))) U = ham_U + ham_dUedge ! set edge site strength
                     if (Hubbard_Mz) then
                        g = im*sqrt(U*Dtau/2)
                     else
                        g = (-1)**(nf + 1)*sqrt(U*Dtau/2)
                     end if
                     Op_V(j, nf)%g = g
                     Op_V(j, nf)%alpha = 0
                     Op_V(j, nf)%type = 2
                     Call Op_set(Op_V(j, nf))
                  end if
               end if
            End do
         end do
      end if
   end subroutine hubbard_V

!--------------------------------------------------------------------
!> Specifiy the equal time and time displaced observables
!--------------------------------------------------------------------
   Subroutine Alloc_obs(Ltau)

      Implicit none
      !>  Ltau=1 if time displaced correlations are considered.
      Integer, Intent(In) :: Ltau
      Integer :: i, j, N, Nt
      integer, allocatable :: ia(:)
      Character(len=64) :: Filename
      character(len=64) :: description(1)
      Character(len=2) :: Channel

      ! Scalar observables
      if (cheatcode(1) > 0) os_last = cheatcode(1)
      Allocate (Obs_scal(os_last))
      Do I = 1, Size(Obs_scal, 1)
         description = ''
         select case (I)
         case (os_kin)
            N = 1; Filename = "kin"  ! Kinetic
         case (os_pot)
            N = 1; Filename = "pot"  ! interaction
         case (os_part)
            N = 1; Filename = "part" ! density
         case (os_ener)
            N = 1; Filename = "ener" ! total
         case (os_hop)
            N = 1; Filename = "hop"  ! nn hopping
         case (os_pair)
            N = 1; Filename = "pair"  ! nn pairing
         case (os_hubd)
            N = 3; Filename = "hubd" ! hubbard
         case (os_deltaS2)
            N = ceiling(L2/2.0); Filename = "deltaS2" ! deltaS - deltaS summed over 2D, of different width
         case (os_gr12)
            N = 16; Filename = "gr12" ! Green's function between cells 1,2
         case (os_deltaS12)
            N = 1; Filename = "deltaS12" ! SS correlation between cells 1,2
         case (os_mstat)
            N = 3; Filename = "mstat" ! auxillary field magnetization: [m, m2, m4, rginv]
         case (os_mstat2)
            N = 8; Filename = "mstat2" ! auxillary field magnetization: [m, m2, ..., m6]
         case (os_mall)
            N = L1*L2*L2; Filename = "mall" ! auxillary field magnetization
         case (os_deltaS1h)
            N = L1*L2*L2; Filename = "DeltaS1h" !
            description = 'DeltaSDeltaS1, assume x periodic'
         case (os_deltaS2h)
            N = L1*L2*L2; Filename = "DeltaS2h" !
            description = 'DeltaSDeltaS2, assume x periodic'
         case (os_deltaS3h)
            N = L1*L2*L2; Filename = "DeltaS3h" !
            description = 'DeltaSDeltaS3, assume x periodic'
         case (os_deltaS4h)
            N = L1*L2*L2; Filename = "DeltaS4h" !
            description = 'DeltaSDeltaS4, assume x periodic'
         case (os_spinxh)
            N = L1*L2*L2; Filename = "SxSxh" !
            description = 'SxSx, assume x periodic'
         case (os_spinyh)
            N = L1*L2*L2; Filename = "SySyh" !
            description = 'SySy, assume x periodic'
         case (os_spinzh)
            N = L1*L2*L2; Filename = "SzSzh" !
            description = 'SzSz, assume x periodic'
         case (os_denh)
            N = L1*L2*L2; Filename = "Denh" !
            description = '(n-1/2)(n-1/2), assume x periodic'
         case (os_grall)
            N = Ndim**2; Filename = "GRall" !
         case (os_deltaSall)
            N = latt%N**2; Filename = "DeltaSall" !
         case default
            Write (output_unit, *) ' Error in Alloc_obs obs_scal', I, os_last, Size(Obs_scal, 1)
         end select
         Call Obser_Vec_make(Obs_scal(I), N, Filename, description=description)
      end do

      ia = unique(mod([integer :: 0, L2/2, L2/4, L2 - 1, L2/2 + 1, L2*3/4], L2))
      interYs = ia(1:min(size(ia), L1))
      open (newunit=j, action='write', file='inter_wires_info')
      write (j, *) 'Recording interwires at Y0='
      write (j, *) interYs
      close (j)
      write (*, *) 'Recording interwires at Y0=', interYs

      ! Equal time correlators
      Allocate (Obs_eq(oe_last))
      Do I = 1, Size(Obs_eq, 1)
         select case (I)
         case (oe_deltas_intra)
            Filename = "deltaS_intrawire"
         case (oe_deltas_inter)
            Filename = "deltaS_interwire"
         case (oe_gr_intra)
            Filename = "gr_intrawire"
         case (oe_gr_inter)
            Filename = "gr_interwire"
         case default
            Write (6, *) "Error in Alloc_obs, obs_eq"
         end select
         Nt = 1
         Channel = "--"
         Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
      end do

      If (Ltau == 1) then
         ! time-displaced correlators
         Allocate (Obs_tau(ot_last))
         Do I = 1, Size(Obs_tau, 1)
            select case (I)
            case (ot_gr_intrawire)
               Channel = 'P'; Filename = "GRT_intrawire"
            case (ot_gr_local_edge)
               Channel = 'P'; Filename = "GRT_localedge"
            case (ot_gr_local_bulk)
               Channel = 'P'; Filename = "GRT_localbulk"
            case default
               Write (6, *) ' Error in Alloc_obs obs_tau'
            end select
            Nt = Ltrot + 1 - 2*Thtrot
            If (Projector) Channel = 'T0'
            Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
         end do
      end if

   contains
      pure function unique(ia) result(ra)
         integer :: i, j, k, l, m
         integer, intent(in) :: ia(:)
         integer, allocatable :: ra(:)
         logical :: wa(size(ia))

         wa = .true.

         do i = 1, size(ia)
            do j = i + 1, size(ia)
               if (i /= j .and. ia(i) == ia(j)) then
                  wa(j) = .false.
               end if
            end do
         end do
         ra = pack(ia, wa)
      end function unique
   End Subroutine Alloc_obs

!--------------------------------------------------------------------
!> @brief
!> Computes equal time observables
!> @details
!> @param [IN] Gr   Complex(:,:,:)
!> \verbatim
!>  Green function: Gr(I,J,nf) = <c_{I,nf } c^{dagger}_{J,nf } > on time slice ntau
!> \endverbatim
!> @param [IN] Phase   Complex
!> \verbatim
!>  Phase
!> \endverbatim
!> @param [IN] Ntau Integer
!> \verbatim
!>  Time slice
!> \endverbatim
!-------------------------------------------------------------------
   subroutine Obser(GR, Phase, Ntau, Mc_step_weight)

      !Use Predefined_Obs

      Implicit none

      Complex(Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim, Ndim, N_FL)
      Complex(Kind=Kind(0.d0)), Intent(IN) :: PHASE
      Integer, INTENT(IN) :: Ntau
      Real(Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight

      !Local
      Complex(Kind=Kind(0.d0)) :: GRC(Ndim, Ndim, N_FL), ZK
      Complex(Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP, ZS, ZZ, ZXY, ZDen, zhop, zpair, zhub
      Integer :: I, J, k, imj, nf, Ix, Iy, nt, ii, jj
      real(rk) :: deltaS14(2, 2)
      Real(Kind=Kind(0.d0)) :: X
      integer :: a, b
      real(8), allocatable, save :: f(:)
      real(8) :: m1, m2, m4

      ZP = PHASE/Real(Phase, kind(0.D0))
      ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
      ZS = ZS*Mc_step_weight

      ! Compute scalar observables.
      Obs_scal(:)%N = Obs_scal(:)%N + 1
      Obs_scal(:)%Ave_sign = Obs_scal(:)%Ave_sign + Real(ZS, kind(0.d0))

      if (obs_scal(1)%N == 1) then
         block
            character(8) :: date
            character(10) :: time
            character(5) :: zone
            call date_and_time(date=date, time=time, zone=zone)
            write (*, *) 'Starting new bin at ', date, ' ', time, ' ', zone
            !call control_print(Group_Comm, '', 'verbose_info')
         end block
      end if

      ! block ! count how often the function is called
      !    integer, save :: ncalled = 0
      !    ncalled = ncalled + 1
      !    if (mod(ncalled, 2000) == 0) then
      !       call toc
      !       write (error_unit, *) 'INFO: ncalled obser', ncalled, Obs_scal(1)%N
      !    end if
      ! end block

      ! initialize GRC
      Do I = 1, Ndim
         Do J = 1, Ndim
            GRC(I, J, :) = -GR(J, I, :)
         End do
         GRC(I, I, :) = 1.D0 + GRC(I, I, :)
      End do
      ! GRC(i,j,nf) = < c^{dagger}_{i,nf } c_{j,nf } >

      ! kinetic energies
      zhop = op_t_energy(hopmat, GRC)/2 ! divide by 2 for nambu
      zpair = op_t_energy(pairmat, GRC)/2 ! divide by 2 for nambu

      Obs_scal(os_hop)%Obs_vec(1) = Obs_scal(os_hop)%Obs_vec(1) + zhop*ZP*ZS
      Obs_scal(os_pair)%Obs_vec(1) = Obs_scal(os_pair)%Obs_vec(1) + zpair*ZP*ZS
      ! kinetic energy
      Obs_scal(os_kin)%Obs_vec(1) = Obs_scal(os_hop)%Obs_vec(1) + Obs_scal(os_pair)%Obs_vec(1)

      ! density
      Zrho = 0
      Do I = 1, Ndim
         if (list(i, 2) <= latt_unit%Norb/2) then ! particle sector
            Zrho = Zrho + sum(Grc(i, i, :))/2
         else ! hole sector
            Zrho = Zrho + sum(Gr(i, i, :))/2
         end if
      end do
      Obs_scal(os_part)%Obs_vec(1) = Obs_scal(os_part)%Obs_vec(1) + Zrho*ZP*ZS

      call hubbard_energy(Latt, Latt_unit, List, GR, GRC, ZS, ZP, Obs_scal(os_hubd))
      Obs_scal(os_pot)%Obs_vec(1) = Obs_scal(os_hubd)%Obs_vec(1)
      ! total energy
      Obs_scal(os_ener)%Obs_vec(1) = Obs_scal(os_kin)%Obs_vec(1) + Obs_scal(os_pot)%Obs_vec(1)

      obs_scal(os_gr12)%obs_vec = obs_scal(os_gr12)%obs_vec + ZP*ZS*reshape(gr([1, 2, 3, 4], [1, 2, 3, 4], 1), [16])

      block
         complex(ck) :: deltas
         deltas = DeltaSDeltaS(GR(1, 3, 1), GR(1, 4, 1), GR(2, 3, 1), GR(2, 4, 1))
         obs_scal(os_deltaS12)%obs_vec(1) = obs_scal(os_deltaS12)%obs_vec(1) + ZP*ZS*deltas
      end block

      call calc_deltaS2

      call all_GR_SS

      ! equal-time Greens functions

      Obs_eq(:)%N = Obs_eq(:)%N + 1
      Obs_eq(:)%Ave_sign = Obs_eq(:)%Ave_sign + Real(ZS, kind(0.d0))

      if (oe_last >= oe_deltas_intra) then
         call correlation_intra_wire(Latt, Latt_unit, List, GR, GRC, ZS, ZP, Obs_eq(oe_deltas_intra))
      end if
      if (oe_last >= oe_deltas_inter) then
         call correlation_inter_wire_pairs(Latt, Latt_unit, List, interYs, GR, GRC, ZS, ZP, Obs_eq(oe_deltas_inter))
      end if
      call gr_intra_wire(Latt, Latt_unit, List, GR, ZS, ZP, 0, Obs_eq(oe_gr_intra))
      call gr_inter_wire(Latt, Latt_unit, List, interYs, GR, ZS, ZP, 0, Obs_eq(oe_gr_inter))

      if (os_last >= os_mstat) then
         f = obser_nsigma(ntau)
         m1 = sum(f)/L1/L2; 
         Obs_scal(os_mstat)%Obs_vec(1) = Obs_scal(os_mstat)%Obs_vec(1) + ZP*ZS*m1
         Obs_scal(os_mstat)%Obs_vec(2) = Obs_scal(os_mstat)%Obs_vec(2) + ZP*ZS*m1**2
         Obs_scal(os_mstat)%Obs_vec(3) = Obs_scal(os_mstat)%Obs_vec(3) + ZP*ZS*m1**4

         Obs_scal(os_mstat2)%Obs_vec = Obs_scal(os_mstat2)%Obs_vec + ZP*ZS*m1**[(i, i=1, 8)]

         if (os_last >= os_mall) then
            do concurrent(i=0:L1 - 1, j=0:L2 - 1, k=0:L2 - 1)
               ii = 1 + i + j*L1 + k*L1*L2
               do jj = 0, L1 - 1
                  a = invlist(latt%invlist(jj, j), 1)
                  b = invlist(latt%invlist(mod(jj + i, L1), k), 1)
                  Obs_scal(os_mall)%Obs_vec(ii) = Obs_scal(os_mall)%Obs_vec(ii) &
                    & + ZP*ZS/L1*f(a)*f(b)
               end do
            end do
         end if
      end if

      do concurrent(i=0:L1 - 1, j=0:L2 - 1, k=0:L2 - 1)
         ii = 1 + i + j*L1 + k*L1*L2
         do jj = 0, L1 - 1
            a = invlist(latt%invlist(jj, j), 1)
            b = invlist(latt%invlist(mod(jj + i, L1), k), 1)

            if (os_last >= os_spinzh) then
               Obs_scal(os_spinxh)%Obs_vec(ii) = Obs_scal(os_spinxh)%Obs_vec(ii) &
               & + ZP*ZS/L1*SxSx(GR(a, b, 1), GR(a + 1, b, 1), GR(a, b + 1, 1), GR(a + 1, b + 1, 1))
               Obs_scal(os_spinyh)%Obs_vec(ii) = Obs_scal(os_spinyh)%Obs_vec(ii) &
               & + ZP*ZS/L1*SySy(GR(a, b, 1), GR(a + 1, b, 1), GR(a, b + 1, 1), GR(a + 1, b + 1, 1))
               Obs_scal(os_spinzh)%Obs_vec(ii) = Obs_scal(os_spinzh)%Obs_vec(ii) &
               & + ZP*ZS/L1*SzSz(GR(a, b, 1), GR(a + 1, b, 1), GR(a, b + 1, 1), GR(a + 1, b + 1, 1))
            end if
            if (os_last >= os_deltaS4h) then
               deltaS14 = DeltaSDeltaS_2x2(GR(a, b, 1), GR(a + 1, b, 1), GR(a, b + 1, 1), GR(a + 1, b + 1, 1))
               Obs_scal(os_deltaS1h)%Obs_vec(ii) = Obs_scal(os_deltaS1h)%Obs_vec(ii) &
                                                   + ZP*ZS/L1*deltaS14(1, 1)
               Obs_scal(os_deltaS2h)%Obs_vec(ii) = Obs_scal(os_deltaS2h)%Obs_vec(ii) &
                                                   + ZP*ZS/L1*deltaS14(2, 1)
               Obs_scal(os_deltaS3h)%Obs_vec(ii) = Obs_scal(os_deltaS3h)%Obs_vec(ii) &
                                                   + ZP*ZS/L1*deltaS14(1, 2)
               Obs_scal(os_deltaS4h)%Obs_vec(ii) = Obs_scal(os_deltaS4h)%Obs_vec(ii) &
                                                   + ZP*ZS/L1*deltaS14(2, 2)
            end if
            if (os_last >= os_denh) then
               Obs_scal(os_denh)%Obs_vec(ii) = Obs_scal(os_denh)%Obs_vec(ii) &
               & + ZP*ZS/L1*densitydensity(GR(a, b, 1), GR(a + 1, b, 1), GR(a, b + 1, 1), GR(a + 1, b + 1, 1) &
                                       & , GR(a, a, 1), GR(a + 1, a + 1, 1), GR(b, b, 1), GR(b + 1, b + 1, 1))
            end if
         end do
      end do

   contains
      ! output entire GR array in nf_calc, sorted
      ! also all to all DeltaS correlators
      subroutine all_GR_SS
         integer :: i, j, k, l, p

         if (os_last >= os_grall) then
            obs_scal(os_grall)%obs_vec = obs_scal(os_grall)%obs_vec + &
                                         ZP*ZS*pack(gr(ordered_orbital_list, ordered_orbital_list, 1), .true.)
         end if

         if (os_last >= os_deltaSall) then
            do concurrent(i=1:ndim/2, j=1:ndim/2)
               k = ordered_site_list(i)
               l = ordered_site_list(j)
               p = i + (j - 1)*ndim/2
               obs_scal(os_deltaSall)%obs_vec(p) = obs_scal(os_deltaSall)%obs_vec(p) + &
                     &    ZP*ZS*DeltaSDeltaS(GR(k, l, 1), GR(k, l + 1, 1), GR(k + 1, l, 1), GR(k + 1, l + 1, 1))
            end do
         end if
      end subroutine

      subroutine calc_deltaS2
         integer :: i, j, k, l, p, dy1, dy2
         real :: mid
         mid = (L2 - 1)/2.0
         if (os_last >= os_deltaS2) then
            do i = 1, Ndim, 2
               do j = 1, Ndim, 2
                  k = list(i, 1)
                  l = list(j, 1)

                  dy1 = nint(mid - abs(modulo(latt%list(k, 2), L2) - mid))
                  dy2 = nint(mid - abs(modulo(latt%list(l, 2), L2) - mid))

                  p = min(dy1, dy2) + 1

                  obs_scal(os_deltaS2)%obs_vec(1:p) = obs_scal(os_deltaS2)%obs_vec(1:p) + &
                        &    ZP*ZS*DeltaSDeltaS(GR(i, j, 1), GR(i, j + 1, 1), GR(i + 1, j, 1), GR(i + 1, j + 1, 1))
               end do
            end do
         end if
      end subroutine calc_deltaS2
   end Subroutine Obser

!--------------------------------------------------------------------
!> @brief
!> Computes time displaced  observables
!> @details
!> @param [IN] NT, Integer
!> \verbatim
!>  Imaginary time
!> \endverbatim
!> @param [IN] GT0, GTT, G00, GTT,  Complex(:,:,:)
!> \verbatim
!>  Green functions:
!>  GT0(I,J,nf) = <T c_{I,nf }(tau) c^{dagger}_{J,nf }(0  )>
!>  G0T(I,J,nf) = <T c_{I,nf }(0  ) c^{dagger}_{J,nf }(tau)>
!>  G00(I,J,nf) = <T c_{I,nf }(0  ) c^{dagger}_{J,nf }(0  )>
!>  GTT(I,J,nf) = <T c_{I,nf }(tau) c^{dagger}_{J,nf }(tau)>
!> \endverbatim
!> @param [IN] Phase   Complex
!> \verbatim
!>  Phase
!> \endverbatim
!-------------------------------------------------------------------
   Subroutine ObserT(NT, GT0, G0T, G00, GTT, PHASE, Mc_step_weight)

      !Use Predefined_Obs

      Implicit none

      Integer, INTENT(IN) :: NT
  Complex(Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim, Ndim, N_FL), G0T(Ndim, Ndim, N_FL), G00(Ndim, Ndim, N_FL), GTT(Ndim, Ndim, N_FL)
      Complex(Kind=Kind(0.d0)), INTENT(IN) :: Phase
      Real(Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight

      !Locals
      Complex(Kind=Kind(0.d0)) :: Z, ZP, ZS, ZZ, ZXY, ZDEN
      Real(Kind=Kind(0.d0)) :: X
      Integer :: IMJ, I, J, I1, J1, no_I, no_J

      ZP = PHASE/Real(Phase, kind(0.D0))
      ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
      ZS = ZS*Mc_step_weight

      If (NT == 0) then
         Obs_tau(:)%N = Obs_tau(:)%N + 1
         Obs_tau(:)%Ave_sign = Obs_tau(:)%Ave_sign + Real(ZS, kind(0.d0))
      End if

      call gr_intra_wire(Latt, Latt_unit, List, GT0, ZS, ZP, Nt, Obs_tau(ot_gr_intrawire))
      if (open_boundary) &
         call gr_local_edge(Latt, Latt_unit, List, GT0, ZS, ZP, Nt, Obs_tau(ot_gr_local_edge))
      call gr_local_bulk(Latt, Latt_unit, List, GT0, ZS, ZP, Nt, Obs_tau(ot_gr_local_bulk))

   end Subroutine OBSERT

   !--------------------------------------------------------------------
   !> @brief
   !> Compute kinetic energy from the sparse matrix Hmat encoding
   !>  Op_T hamiltonian
   !--------------------------------------------------------------------
   pure function op_t_energy(Hmat, GRC) result(E)
      real(rk) :: E
      complex(ck) :: tE
      type(Spmatz), intent(in) :: Hmat(:)
      complex(ck), intent(in) :: GRC(:, :, :)
      integer :: i, n

      tE = 0
      do n = 1, size(Hmat, 1)
         do i = 1, (Hmat(n)%nnz)
            tE = tE + GRC(Hmat(n)%entry(i)%i, Hmat(n)%entry(i)%j, n)*Hmat(n)%entry(i)%val
         end do
      end do
      E = real(tE)
   end function op_t_energy

   Subroutine hubbard_energy(Latt, Latt_unit, List, GR, GRC, ZS, ZP, ObsE)

      Type(Lattice), Intent(in) :: Latt
      Type(Unit_cell), Intent(in) :: Latt_unit
      Integer, Intent(In) :: LIST(:, :)

      Complex(Kind=Kind(0.d0)), Intent(In) :: GR(:, :, :), GRC(:, :, :), ZS, ZP
      Type(Obser_Vec), Intent(inout) :: ObsE ! ObsE gives the energy

      ! Local
      integer :: lower_right_uc, upper_right_uc, i, no, k, k1, nh
      real(rk) :: U
      complex(kind=kind(0.d0)) :: hubbardEnergy, hubbardEnergy1, hubbardEnergy2

      hubbardEnergy = 0
      hubbardEnergy1 = 0
      hubbardEnergy2 = 0

      nh = latt_unit%norb/2 ! half of nambu basis length
      ! meaure throughout lattice sites
      ! alternatively one may also scan through Op_V
      do I = 1, Latt%N
         do no = 1, nh
            k1 = invlist(I, no)
            k = k1 + nh
            if (.not. is_excluded_site(k1)) then
               if (is_edge_site(k1)) then
                  U = ham_U + Ham_dUedge
               else
                  U = ham_U
               end if
               hubbardEnergy = hubbardEnergy + U*(GRC(k1, k1, 1) - GRC(k, k, 1))/2*(GRC(k1, k1, 2) - GRC(k, k, 2))/2
            hubbardEnergy1 = hubbardEnergy1 + U*(GR(k, k, 1)*conjg(GRC(k, k, 1)) - 0.5*GR(k, k, 1) - 0.5*conjg(GRC(k, k, 1)) + 0.25)
               hubbardEnergy2 = hubbardEnergy2 + U*(GRC(k1, k1, 1) - 0.5)*(GRC(k1, k1, 2) - 0.5)
            end if
         end do
      end do

      ObsE%Obs_vec(1) = ObsE%Obs_vec(1) + ZP*ZS*hubbardEnergy
      ObsE%Obs_vec(2) = ObsE%Obs_vec(2) + ZP*ZS*hubbardEnergy1
      ObsE%Obs_vec(3) = ObsE%Obs_vec(3) + ZP*ZS*hubbardEnergy2

   end Subroutine hubbard_energy

   pure function DeltaSDeltaS(GRmn, GRm2n, GRmn2, GRm2n2)
      complex(ck) :: DeltaSDeltaS
      complex(ck), intent(in) :: GRmn, GRm2n, GRmn2, GRm2n2
      DeltaSDeltaS = conjg(GRmn)*GRmn + conjg(GRm2n)*GRm2n + conjg(GRmn2)*GRmn2 + conjg(GRm2n2)*GRm2n2
   end function DeltaSDeltaS

   pure function DeltaSDeltaS_2x2(GRmn, GRm2n, GRmn2, GRm2n2)
      real(rk) :: DeltaSDeltaS_2x2(2, 2)
      complex(ck), intent(in) :: GRmn, GRm2n, GRmn2, GRm2n2
      DeltaSDeltaS_2x2(1, 1) = real(conjg(GRmn)*GRmn + conjg(GRm2n)*GRm2n + conjg(GRmn2)*GRmn2 + conjg(GRm2n2)*GRm2n2)
      DeltaSDeltaS_2x2(1, 2) = real(conjg(GRmn)*GRmn - conjg(GRm2n)*GRm2n - conjg(GRmn2)*GRmn2 + conjg(GRm2n2)*GRm2n2)
      DeltaSDeltaS_2x2(2, 1) = real(conjg(GRmn)*GRmn - conjg(GRm2n)*GRm2n + conjg(GRmn2)*GRmn2 - conjg(GRm2n2)*GRm2n2)
      DeltaSDeltaS_2x2(2, 2) = real(conjg(GRmn)*GRmn + conjg(GRm2n)*GRm2n - conjg(GRmn2)*GRmn2 - conjg(GRm2n2)*GRm2n2)
   end function DeltaSDeltaS_2x2

   ! SxSx, good for equal time and time displaced
   pure real(rk) function SxSx(GRmn, GRm2n, GRmn2, GRm2n2)
      complex(ck), intent(in) :: GRmn, GRm2n, GRmn2, GRm2n2
      SxSx = real(GRmn*conjg(GRm2n2) - GRm2n*conjg(GRmn2))/2
   end function SxSx

   ! SySy, good for equal time and time displaced
   pure real(rk) function SySy(GRmn, GRm2n, GRmn2, GRm2n2)
      complex(ck), intent(in) :: GRmn, GRm2n, GRmn2, GRm2n2
      SySy = real(GRmn*conjg(GRm2n2) + GRm2n*conjg(GRmn2))/2
   end function SySy

   ! SzSz, good for equal time and time displaced
   pure real(ck) function SzSz(GRmn, GRm2n, GRmn2, GRm2n2)
      complex(ck), intent(in) :: GRmn, GRm2n, GRmn2, GRm2n2
      SzSz = real(GRmn*GRm2n2)/2
   end function SzSz

   ! density-density, good for equal time *only*
   pure complex(ck) function densitydensity(GRmn, GRm2n, GRmn2, GRm2n2, GRmm, GRm2m2, GRnn, GRn2n2)
      complex(ck), intent(in) :: GRmn, GRm2n, GRmn2, GRm2n2
      complex(ck), intent(in) :: GRmm, GRm2m2, GRnn, GRn2n2
      densitydensity = GRmn*GRm2n2 - GRmn2*GRm2n + (GRmm - GRm2m2)*(GRnn - GRn2n2)/4
   end function densitydensity

   Subroutine correlation_intra_wire(Latt, Latt_unit, List, GR, GRC, ZS, ZP, Obs)
      Type(Lattice), Intent(in) :: Latt
      Type(Unit_cell), Intent(in) :: Latt_unit
      Integer, Intent(In) :: LIST(:, :)
      !Real (Kind=Kind(0.d0)), Intent(In)    :: Y0
      Complex(Kind=Kind(0.d0)), Intent(In) :: GR(:, :, :), GRC(:, :, :), ZS, ZP
      Type(Obser_Latt), Intent(inout) :: Obs   ! Obs(Latt%L2), each gives the correlation function, ObsX%Obs_Latt(imj,tau,Ix,Iy)

      ! Local
      Integer :: X1, X2, I1, I2, a1, a2, m, n, imj, LY    ! LY is the position of the layer, Latt%L2 is the number of layer
      real(kind=kind(0.0d0)) :: s

      ! meaure, only 1st layer / one of the nambu
      do LY = 0, L2 - 1
         do X1 = 1, L1
            I1 = latt%invlist(X1, LY)
            do X2 = 1, L1
               I2 = latt%invlist(X2, LY)
               imj = latt%imj(I1, latt%invlist(X2, 0))

               m = invlist(I1, 1)
               n = invlist(I2, 1)
               ! note that m,n are the particles, and m+1,n+1 are the holes or nambu partners
               Obs%Obs_Latt(imj, 1, :, :) = Obs%Obs_Latt(imj, 1, :, :) + &
                   &  ZP*ZS*DeltaSDeltaS_2x2(GR(m, n, 1), GR(m, n + 1, 1), GR(m + 1, n, 1), GR(m + 1, n + 1, 1))*L2

            end do
         end do
      end do
   end Subroutine correlation_intra_wire

   Subroutine correlation_inter_wire_pairs(Latt, Latt_unit, List, Y0s, GR, GRC, ZS, ZP, obs)

      Type(Lattice), Intent(in) :: Latt
      Type(Unit_cell), Intent(in) :: Latt_unit
      Integer, Intent(In) :: LIST(:, :)
      Integer, Intent(In) :: Y0s(:)          ! Y0 is the position of edge layer
      Complex(Kind=Kind(0.d0)), Intent(In) :: GR(:, :, :), GRC(:, :, :), ZS, ZP
      Type(Obser_latt), Intent(inout) :: Obs ! obs(Latt%L2_p), each gives the layer susceptibility function from boundary to any layers

      ! Local
      Integer :: N_FL, X0, dY, Y0, X1, Y1, a1, a2, I1, I2, m, n, imj, iY
      real(rk) :: s

      do iY = 1, size(Y0s)
         Y0 = Y0s(iY)
         do dY = 0, L2 - 1 ! distance between two layers
            Y1 = mod(dY + Y0, L2)
            imj = latt%imj(latt%invlist(iY - 1, dY), latt%invlist(0, 0))

            do X0 = 1, L1        ! the unit cell at the edge
               I1 = Latt%invlist(X0, 0)
               do X1 = 1, L1         ! unit cell at the layer YL
                  I2 = Latt%invlist(X1, Y1)

                  m = invlist(I1, 1)
                  n = invlist(I2, 1)

                  Obs%Obs_Latt(imj, 1, :, :) = Obs%Obs_Latt(imj, 1, :, :) + &
                      &  ZP*ZS*DeltaSDeltaS_2x2(GR(m, n, 1), GR(m, n + 1, 1), GR(m + 1, n, 1), GR(m + 1, n + 1, 1))/L1/L1
               end do
            end do
         end do
      end do

   end Subroutine correlation_inter_wire_pairs

   Subroutine gr_intra_wire(Latt, Latt_unit, List, GT0, ZS, ZP, Nt, Obs)
      Type(Lattice), Intent(in) :: Latt
      Type(Unit_cell), Intent(in) :: Latt_unit
      Integer, Intent(In) :: LIST(:, :)
      !Real (Kind=Kind(0.d0)), Intent(In)    :: Y0
      Integer, INTENT(IN) :: NT
      Complex(Kind=Kind(0.d0)), INTENT(IN) :: GT0(:, :, :)
      Complex(Kind=Kind(0.d0)), Intent(In) :: ZS, ZP
      Type(Obser_Latt), Intent(inout) :: Obs   ! Obs(Latt%L2), each gives the correlation function, ObsX%Obs_Latt(imj,tau,Ix,Iy)

      ! Local
      Integer :: X1, X2, I1, I2, a1, a2, m, n, imj, LY    ! LY is the position of the layer, Latt%L2 is the number of layer
      real(kind=kind(0.0d0)) :: s

      ! meaure, only 1st layer / one of the nambu
      do LY = 0, L2 - 1
         do X1 = 1, L1
            I1 = latt%invlist(X1, LY)
            do X2 = 1, L1
               I2 = latt%invlist(X2, LY)
               imj = latt%imj(I1, latt%invlist(X2, 0))
               do a1 = 1, latt_unit%norb ! to denote sublattice, not Nambu, Nambu is summed in the following
                  do a2 = 1, latt_unit%norb
                     m = invlist(I1, a1)
                     n = invlist(I2, a2)
                     ! note that m,n are the particles, and m+1,n+1 are the holes or nambu partners
                     Obs%Obs_Latt(imj, Nt + 1, a1, a2) = Obs%Obs_Latt(imj, Nt + 1, a1, a2) + &
                         &  ZP*ZS*L2*GT0(m, n, 1)
                  end do
               end do
            end do
         end do
      end do
   end Subroutine gr_intra_wire

   Subroutine gr_inter_wire(Latt, Latt_unit, List, Y0s, GT0, ZS, ZP, Nt, Obs)
      Type(Lattice), Intent(in) :: Latt
      Type(Unit_cell), Intent(in) :: Latt_unit
      Integer, Intent(In) :: LIST(:, :)
      Integer, Intent(In) :: Y0s(:)          ! Y0 is the position of edge layer
      integer, intent(in) :: Nt
      Complex(Kind=Kind(0.d0)), Intent(In) :: GT0(:, :, :), ZS, ZP
      Type(Obser_latt), Intent(inout) :: Obs

      ! Local
      Integer :: N_FL, X0, dY, Y0, X1, Y1, a1, a2, I1, I2, m, n, imj, iY
      real(rk) :: s

      do iY = 1, size(Y0s)
         Y0 = Y0s(iY)
         do dY = 0, L2 - 1 ! distance between two layers
            Y1 = mod(dY + Y0, L2)
            imj = latt%imj(latt%invlist(iY - 1, dY), latt%invlist(0, 0)) ! iY'th index at dX, dY index at dY
            do X0 = 1, L1        ! the unit cell at the edge
               I1 = Latt%invlist(X0, 0)
               do X1 = 1, L1         ! unit cell at the layer YL
                  I2 = Latt%invlist(X1, Y1)
                  do a1 = 1, latt_unit%Norb         ! to denote sublattice, not Nambu, Nambu is summed in the following
                     do a2 = 1, latt_unit%Norb
                        m = invlist(I1, a1)
                        n = invlist(I2, a2)
                        Obs%Obs_Latt(imj, Nt + 1, a1, a2) = Obs%Obs_Latt(imj, Nt + 1, a1, a2) + &
                                 &  ZP*ZS*GT0(m, n, 1)*L2/L1
                     end do
                  end do
               end do
            end do
         end do
      end do
   end Subroutine gr_inter_wire

   Subroutine gr_local_edge(Latt, Latt_unit, List, GT0, ZS, ZP, Nt, Obs)
      Type(Lattice), Intent(in) :: Latt
      Type(Unit_cell), Intent(in) :: Latt_unit
      Integer, Intent(In) :: LIST(:, :)
      !Real (Kind=Kind(0.d0)), Intent(In)    :: Y0
      Integer, INTENT(IN) :: NT
      Complex(Kind=Kind(0.d0)), INTENT(IN) :: GT0(:, :, :)
      Complex(Kind=Kind(0.d0)), Intent(In) :: ZS, ZP
      Type(Obser_Latt), Intent(inout) :: Obs   ! Obs(Latt%L2), each gives the correlation function, ObsX%Obs_Latt(imj,tau,Ix,Iy)

      ! Local
      Integer :: I1, I2
      integer, allocatable, save :: edge_sites(:) ! sites belonging to the edge cells, particles only

      if (.not. allocated(edge_sites)) then
         edge_sites = invlist(edge_cell_list, 1)
      end if
      do i1 = 1, size(edge_sites)
         i2 = edge_sites(i1)
         Obs%Obs_Latt(:, Nt + 1, :, :) = Obs%Obs_Latt(:, Nt + 1, :, :) + ZP*ZS*GT0(i2, i2, 1)/size(edge_sites)
      end do
   end Subroutine gr_local_edge

   Subroutine gr_local_bulk(Latt, Latt_unit, List, GT0, ZS, ZP, Nt, Obs)
      Type(Lattice), Intent(in) :: Latt
      Type(Unit_cell), Intent(in) :: Latt_unit
      Integer, Intent(In) :: LIST(:, :)
      !Real (Kind=Kind(0.d0)), Intent(In)    :: Y0
      Integer, INTENT(IN) :: NT
      Complex(Kind=Kind(0.d0)), INTENT(IN) :: GT0(:, :, :)
      Complex(Kind=Kind(0.d0)), Intent(In) :: ZS, ZP
      Type(Obser_Latt), Intent(inout) :: Obs   ! Obs(Latt%L2), each gives the correlation function, ObsX%Obs_Latt(imj,tau,Ix,Iy)

      ! Local
      Integer :: I1, I2
      integer, allocatable, save :: center_sites(:) ! sites at the center

      if (.not. allocated(center_sites)) then
         if (mod(L2, 2) == 0) then ! even L2, center is between i2=L2/2-1'th and L2/2'th (note integer division)
            allocate (center_sites(2*L1))
            center_sites(1:L1) = pack(invlist(latt%invlist([(i1, i1=1, L1)], L2/2 - 1), 1), .true.)
            center_sites((L1 + 1):(2*L1)) = pack(invlist(latt%invlist([(i1, i1=1, L1)], L2/2), 1), .true.)
         else ! odd L2, center is at i2=L2/2'th row
            center_sites = pack(invlist(latt%invlist([(i1, i1=1, L1)], L2/2), 1), .true.)
         end if
      end if

      do i1 = 1, size(center_sites)
         i2 = center_sites(i1)
         Obs%Obs_Latt(:, Nt + 1, :, :) = Obs%Obs_Latt(:, Nt + 1, :, :) + ZP*ZS*GT0(i2, i2, 1)/size(center_sites)
      end do
   end Subroutine gr_local_bulk

   ! function to return current values of auxiliary fields
   function obser_nsigma(ntau) result(f)
      integer, intent(in) :: ntau
      real(8), allocatable :: f(:)
      real(8), parameter :: x1 = sqrt(2*(3 - sqrt(6d0))), x2 = sqrt(2*(3 + sqrt(6d0)))

      f = sign(merge(x1, x2, abs(nsigma%f(:, ntau)) == 1), nsigma%f(:, ntau))
   end function

!--------------------------------------------------------------------
!> @brief
!> Reconstructs dependent flavors of the configuration's weight.
!> @details
!> This has to be overloaded in the Hamiltonian submodule.
!--------------------------------------------------------------------
   subroutine weight_reconstruction(weight)
      implicit none
      complex(Kind=Kind(0.d0)), Intent(inout) :: weight(:)

      weight(nf_reconst) = 1
      weight(nf_calc) = abs(Weight(nf_calc)) ! nambu
      !weight(nf_reconst) = abs(Weight(nf_calc))**2/Weight(nf_calc) ! no nambu

   end subroutine weight_reconstruction

!--------------------------------------------------------------------
!> @brief
!> Reconstructs dependent flavors of equal time Greens function
!> @details
!> This has to be overloaded in the Hamiltonian submodule.
!> @param [INOUT] Gr   Complex(:,:,:)
!> \verbatim
!>  Green function: Gr(I,J,nf) = <c_{I,nf } c^{dagger}_{J,nf } > on time slice ntau
!> \endverbatim
!-------------------------------------------------------------------
   subroutine GR_reconstruction(GR)

      Implicit none

      Complex(Kind=Kind(0.d0)), INTENT(INOUT) :: GR(Ndim, Ndim, N_FL)
      Integer :: I, J, imj
      real(kind=kind(0.d0)) :: X, ZZ

      if (N_Fl == 1) return
      Do J = 1, Ndim
         Do I = 1, Ndim
            GR(I, J, nf_reconst) = conjg(GR(I, J, nf_calc))
         End do
      End do

   end Subroutine GR_reconstruction

!--------------------------------------------------------------------
!> @brief
!> Reconstructs dependent flavors of time displaced Greens function G0T and GT0
!> @details
!> This has to be overloaded in the Hamiltonian submodule.
!> @param [INOUT] GT0, G0T,  Complex(:,:,:)
!> \verbatim
!>  Green functions:
!>  GT0(I,J,nf) = <T c_{I,nf }(tau) c^{dagger}_{J,nf }(0  )>
!>  G0T(I,J,nf) = <T c_{I,nf }(0  ) c^{dagger}_{J,nf }(tau)>
!> \endverbatim
!-------------------------------------------------------------------
   Subroutine GRT_reconstruction(GT0, G0T)
      Implicit none

      Complex(Kind=Kind(0.d0)), INTENT(INOUT) :: GT0(Ndim, Ndim, N_FL), G0T(Ndim, Ndim, N_FL)
      Integer :: I, J, imj
      real(kind=kind(0.d0)) :: X, ZZ

      if (N_Fl == 1) return
      Do J = 1, Latt%N
         Do I = 1, Latt%N
            G0T(I, J, nf_reconst) = conjg(G0T(I, J, nf_calc))
            GT0(I, J, nf_reconst) = conjg(GT0(I, J, nf_calc))
         end do
      end do

   end Subroutine GRT_reconstruction

   include 'yge5_util.f90'

end submodule ham_TSC_smod
