subroutine Build_S(S,Nr_atoms,HDIM,Max_Nr_Neigh,Element_Type,RX,RY,RZ,LBox,H_INDEX_START,H_INDEX_END, & 
 nrnnlist,nnRx,nnRy,nnRz,nnType)

use omp_lib
implicit none
integer,    parameter      :: PREC = 8
integer,    intent(in)     :: Nr_atoms, HDIM, Max_Nr_Neigh
real(PREC), intent(in)     :: RX(Nr_atoms), RY(Nr_atoms), RZ(Nr_atoms), LBox(3)
integer,    intent(in)     :: H_INDEX_START(Nr_atoms), H_INDEX_END(Nr_atoms)
real(PREC), intent(out)    :: S(HDIM,HDIM)
integer                    :: IDim, JDim 
real(PREC)                 :: Ra(3), Rb(3), diagonal(2), s_0(4,4)
real(PREC)                 :: fss_sigma(14),fsp_sigma(14),fps_sigma(14),fpp_sigma(14),fpp_pi(14)
character(10), intent(in)  :: Element_Type(Nr_atoms)
character(10)              ::  Type_pair(2)

integer,    intent(in)     :: nrnnlist(Nr_atoms), nnType(Nr_atoms,Max_Nr_Neigh)
real(PREC), intent(in)     :: nnRx(Nr_atoms,Max_Nr_Neigh)
real(PREC), intent(in)     :: nnRy(Nr_atoms,Max_Nr_Neigh), nnRz(Nr_atoms,Max_Nr_Neigh)


integer                    :: I,J,II,JJ,II_S,JJ_S, nthreads,threads,id

!write(*,*) ' Build_S RX(1) = ', RX(1)
!nthreads = 8
!call omp_set_num_threads(nthreads)

S = 0.0D0  ! Replace by sparse ELLPACK matrix
s_0 = 0.D0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,II,JJ,II_S,JJ_S,Ra,Rb,fss_sigma,fsp_sigma,fps_sigma,fpp_sigma,fpp_pi) &
!$OMP PRIVATE(Type_pair,diagonal,IDim,JDim,s_0)
do I = 1,Nr_atoms
!    threads = omp_get_num_threads()
!    id = omp_get_thread_num()
!    write(*,*) 'NUM THREADS:', threads,id, omp_get_max_threads()
  Type_pair(1) = Element_Type(I)
  Ra(1) = RX(I)
  Ra(2) = RY(I)
  Ra(3) = RZ(I)
  IDim = H_INDEX_END(I)-H_INDEX_START(I)+1
  do J = 1,nrnnlist(I)
    Type_pair(2) = Element_Type(nnType(I,J))         
    Rb(1) = nnRx(I,J)
    Rb(2) = nnRy(I,J)
    Rb(3) = nnRz(I,J)
    JDim = H_INDEX_END(nnType(I,J))-H_INDEX_START(nnType(I,J))+1
! Overlap block for a-b atom pair
    call LoadBondIntegralParameters_S(fss_sigma,fsp_sigma,fps_sigma,fpp_sigma,fpp_pi,Type_pair)
    diagonal(1:2) = [1,1]
    !call Slater_Koster_Block(s_0,Ra,Rb,LBox,Type_pair,fss_sigma,fsp_sigma,fps_sigma,fpp_sigma,fpp_pi,diagonal)
    call Slater_Koster_Pair(s_0,Ra,Rb,LBox,Type_pair,fss_sigma,fsp_sigma,fps_sigma,fpp_sigma,fpp_pi,diagonal)
    do II = 1,IDim
      II_S = H_INDEX_START(I) + II - 1
      do JJ = 1,JDim
        JJ_S = H_INDEX_START(nnType(I,J)) + JJ - 1
        S(II_S,JJ_S) = s_0(II,JJ) ! Replace by sparse ELLPACK matrix
        S(JJ_S,II_S) = s_0(II,JJ)
      enddo
    enddo

  enddo
enddo
!$OMP END PARALLEL DO

end subroutine Build_S
