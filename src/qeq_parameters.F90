module qeq_parameters
  
  IMPLICIT NONE

  ! ev = 1.602e-19 J
  ! K = 1.38xe-23 J/K
  ! mass = 1.67262192e-27 Kg
  ! E_ket = mv^2 = kg * (m/s)^2 = 1.672e-27 au * (ang/fs)^2 e10 J
  !              = 1.672au * (ang/fs) e-17/1.602e-19 eV
  !              = 167.xx / 1.602xx  eV
  !
  ! Force = kg*m/s = 1/1.6605e-27 au * 1.0e-5 ang/fs
  !       = J/m = 1.602e-19 ev/e10 ang 
  !       1.602e-29 ev/ang, F/m -> 1.602/1.6605 e-2, (ev/ang/au)

  ! mv^2 = J = kg * (m/s)^2 = kg * (ang/fs) * e-20
  ! mv^2 = kg * (ang/fs) * e-20 / 1.602e-19 

  integer, parameter       :: PREC = 8
  real(PREC), parameter    :: F2V = 0.01602176487D0/1.660548782D0  !
  real(PREC), parameter    :: MVV2KE = 166.0538782D0/1.602176487D0 !
  real(PREC), parameter    :: KE2T = 1.D0/0.000086173435D0         ! ev to Kelvin
  CHARACTER(LEN = 1000) :: OUTFILE = "log.latte"
  real(PREC) :: qeq_xi(300)

  real(PREC) :: Coeff(7), C0, C1, C2, C3, C4, C5, C6, kappa, Time, EKIN, Energy
  integer :: HDIM, Nocc

  character(1) :: PairPotType1(10), PairPotType2(10)
  real(PREC) :: LBox(3), COULCUT
  real(PREC) :: SCF_ERR
  real(PREC) :: Coulomb_Force_Real_I(3), Coulomb_Pot_Real_I
  real(PREC) :: proj_tmp, proj_tmp2
  real(PREC) :: First, lambda, temp2, temp, RMSE, Fel, Temperature
  real(PREC) :: ERep, ECoul, EBand, EEnt, EPOT

logical :: exact_solution, LATTEINEXISTS
logical :: printcharges
logical :: printpot

character(10), allocatable :: ATELE(:)
integer, allocatable :: H_INDEX_START(:), H_INDEX_END(:), NrOrb(:)
integer, allocatable :: nrnnlist(:), nnType(:,:)
integer, allocatable :: nnStruct(:,:), nrnnStruct(:)
real(PREC), allocatable :: RX(:), RY(:), RZ(:), q(:), qx(:), qqx(:)
real(PREC), allocatable :: Znuc(:), Mnuc(:), Hubbard_U(:), Xi(:), bb(:)
real(PREC), allocatable :: CC(:,:), AA(:,:), xx(:)
real(PREC), allocatable :: AAI(:,:), KK(:,:)
real(PREC), allocatable :: nndist(:,:), nnRx(:,:)
real(PREC), allocatable :: nnRy(:,:), nnRz(:,:), Coulomb_Pot_k(:)
real(PREC), allocatable :: vi(:,:), v(:), q_tmp(:), dq_dv(:)
real(PREC), allocatable :: fi(:,:), vK(:),vJ(:)
real(PREC), allocatable :: dr(:), IdentRes(:), II(:,:)
real(PREC), allocatable :: OO(:,:),MM(:,:)
real(PREC), allocatable :: Coulomb_Pot_Real(:), Coulomb_Force_Real(:,:)
real(PREC), allocatable :: Coulomb_Force_k(:,:), Coulomb_Pot(:), Coulomb_Force(:,:)
real(PREC), allocatable :: JJ(:,:), KK0(:,:), DD(:,:)
real(PREC), allocatable :: gq(:), KRes(:), S_Ent, VV(:), PairForces(:,:)
real(PREC), allocatable :: SKForce(:,:), FPUL(:,:), FSCOUL(:,:), FTOT(:,:)
real(PREC), allocatable :: n(:),n_0(:),n_1(:),n_2(:),n_3(:),n_4(:)
real(PREC), allocatable :: n_5(:), n_6(:)
real(PREC), allocatable :: VX(:), VY(:), VZ(:)
real(PREC), allocatable :: dn2dt2(:), Res(:), K0Res(:)

real(PREC) :: alpha, beta, mu0, mu1
real(PREC), allocatable :: D_atomic(:), H0(:,:), H(:,:), H1(:,:)
real(PREC), allocatable :: D0(:,:), D(:,:)
real(PREC), allocatable :: S(:,:), Z(:,:), X(:,:), Y(:,:) 
real(PREC), allocatable :: HOrth(:,:)
real(PREC), allocatable :: dSx(:,:), dSy(:,:), dSz(:,:)
real(PREC), allocatable :: dHx(:,:), dHy(:,:), dHz(:,:)
real(PREC), allocatable :: QQ(:,:), ee(:), Fe_vec(:)  ! Eigenvectors, eigenvalues and Fermi occupation
real(PREC), allocatable :: Coulomb_Pot_Exact(:)       !! ANDERS NEW FEB 9 2023

real(PREC) :: PotCoef(16,10)

integer :: VERBOSE

!  qeq_params = 0.D0

  contains

    subroutine get_qeq_params()
    implicit none

    qeq_xi = 0.01D0

    qeq_xi(1) = 0.13D0   ! H electron negativity 
    qeq_xi(16) = -3.5D0  ! O
    qeq_xi(238) = 0.6    ! U

    end subroutine
    
    

end module qeq_parameters
