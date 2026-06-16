program embedding_fortran
  use, intrinsic :: iso_c_binding
  implicit none

  interface
    function otfp_engine_create(ncenters, ncvs, nrestr, seed) bind(C, name="otfp_engine_create")
      import :: c_int, c_long, c_ptr
      integer(c_int), value :: ncenters, ncvs, nrestr
      integer(c_long), value :: seed
      type(c_ptr) :: otfp_engine_create
    end function

    subroutine otfp_engine_free(engine) bind(C, name="otfp_engine_free")
      import :: c_ptr
      type(c_ptr), value :: engine
    end subroutine

    function otfp_engine_add_center(engine, natoms, indices, masses) bind(C, name="otfp_engine_add_center")
      import :: c_int, c_double, c_ptr
      type(c_ptr), value :: engine
      integer(c_int), value :: natoms
      integer(c_int), dimension(*) :: indices
      real(c_double), dimension(*) :: masses
      integer(c_int) :: otfp_engine_add_center
    end function

    function otfp_engine_set_center_position(engine, center, xyz) bind(C, name="otfp_engine_set_center_position")
      import :: c_int, c_double, c_ptr
      type(c_ptr), value :: engine
      integer(c_int), value :: center
      real(c_double), dimension(*) :: xyz
      integer(c_int) :: otfp_engine_set_center_position
    end function

    function otfp_engine_center_force(engine, center) bind(C, name="otfp_engine_center_force")
      import :: c_int, c_ptr
      type(c_ptr), value :: engine
      integer(c_int), value :: center
      type(c_ptr) :: otfp_engine_center_force
    end function

    function otfp_engine_add_cv(engine, cvtype, ncenters, centers, zmin, zmax, boundary, boundary_k) &
        bind(C, name="otfp_engine_add_cv")
      import :: c_char, c_int, c_double, c_ptr
      type(c_ptr), value :: engine
      character(kind=c_char), dimension(*) :: cvtype
      integer(c_int), value :: ncenters
      integer(c_int), dimension(*) :: centers
      real(c_double), value :: zmin, zmax
      character(kind=c_char), dimension(*) :: boundary
      real(c_double), value :: boundary_k
      integer(c_int) :: otfp_engine_add_cv
    end function

    function otfp_engine_add_restraint(engine, k, target, ncvs, coefficients, function_name, &
        zmin, zmax, boundary, boundary_k) bind(C, name="otfp_engine_add_restraint")
      import :: c_char, c_int, c_double, c_ptr
      type(c_ptr), value :: engine
      real(c_double), value :: k, target
      integer(c_int), value :: ncvs
      real(c_double), dimension(*) :: coefficients
      character(kind=c_char), dimension(*) :: function_name
      real(c_double), value :: zmin, zmax
      character(kind=c_char), dimension(*) :: boundary
      real(c_double), value :: boundary_k
      integer(c_int) :: otfp_engine_add_restraint
    end function

    function otfp_engine_compute(engine, first_step, step, external_bias) bind(C, name="otfp_engine_compute")
      import :: c_int, c_double, c_ptr
      type(c_ptr), value :: engine
      integer(c_int), value :: first_step, step
      real(c_double), value :: external_bias
      integer(c_int) :: otfp_engine_compute
    end function

    function otfp_engine_bias_energy(engine) bind(C, name="otfp_engine_bias_energy")
      import :: c_double, c_ptr
      type(c_ptr), value :: engine
      real(c_double) :: otfp_engine_bias_energy
    end function
  end interface

  type(c_ptr) :: engine
  type(c_ptr) :: force_ptr
  real(c_double), pointer :: force(:)
  integer(c_int), target :: ind0(1) = [0_c_int]
  integer(c_int), target :: ind1(1) = [1_c_int]
  integer(c_int), target :: cv_centers(2) = [0_c_int, 1_c_int]
  real(c_double), target :: mass(1) = [1.0_c_double]
  real(c_double), target :: cvc(1) = [1.0_c_double]
  real(c_double), target :: p0(3) = [0.0_c_double, 0.0_c_double, 0.0_c_double]
  real(c_double), target :: p1(3) = [1.0_c_double, 0.0_c_double, 0.0_c_double]
  character(kind=c_char), target :: bond(5) = [character(kind=c_char) :: 'B', 'O', 'N', 'D', c_null_char]
  character(kind=c_char), target :: nada(5) = [character(kind=c_char) :: 'N', 'A', 'D', 'A', c_null_char]
  character(kind=c_char), target :: harmonic(9) = [character(kind=c_char) :: &
      'H', 'A', 'R', 'M', 'O', 'N', 'I', 'C', c_null_char]
  integer(c_int) :: code
  real(c_double) :: energy

  engine = otfp_engine_create(2_c_int, 1_c_int, 1_c_int, 1234_c_long)
  if (.not. c_associated(engine)) error stop 1

  if (otfp_engine_add_center(engine, 1_c_int, ind0, mass) /= 0_c_int) error stop 2
  if (otfp_engine_add_center(engine, 1_c_int, ind1, mass) /= 1_c_int) error stop 3
  if (otfp_engine_add_cv(engine, bond, 2_c_int, cv_centers, 0.0_c_double, 0.0_c_double, &
      nada, 0.0_c_double) /= 0_c_int) error stop 4
  if (otfp_engine_add_restraint(engine, 1.0_c_double, 0.0_c_double, 1_c_int, cvc, harmonic, &
      0.0_c_double, 0.0_c_double, nada, 0.0_c_double) /= 0_c_int) error stop 5

  if (otfp_engine_set_center_position(engine, 0_c_int, p0) /= 0_c_int) error stop 6
  if (otfp_engine_set_center_position(engine, 1_c_int, p1) /= 0_c_int) error stop 7

  code = otfp_engine_compute(engine, 0_c_int, 1_c_int, 0.0_c_double)
  if (code /= 0_c_int) error stop 8

  force_ptr = otfp_engine_center_force(engine, 0_c_int)
  if (.not. c_associated(force_ptr)) error stop 9
  call c_f_pointer(force_ptr, force, [3])

  energy = otfp_engine_bias_energy(engine)
  if (energy /= energy) error stop 10
  if (any(force /= force)) error stop 11

  call otfp_engine_free(engine)
  call otfp_engine_free(c_null_ptr)
end program
