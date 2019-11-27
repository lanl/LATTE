!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 2010.  Los Alamos National Security, LLC. This material was    !
! produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos !
! National Laboratory (LANL), which is operated by Los Alamos National     !
! Security, LLC for the U.S. Department of Energy. The U.S. Government has !
! rights to use, reproduce, and distribute this software.  NEITHER THE     !
! GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY,     !
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS         !
! SOFTWARE.  If software is modified to produce derivative works, such     !
! modified software should be clearly marked, so as not to confuse it      !
! with the version available from LANL.                                    !
!                                                                          !
! Additionally, this program is free software; you can redistribute it     !
! and/or modify it under the terms of the GNU General Public License as    !
! published by the Free Software Foundation; version 2.0 of the License.   !
! Accordingly, this program is distributed in the hope that it will be     !
! useful, but WITHOUT ANY WARRANTY; without even the implied warranty of   !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General !
! Public License for more details.                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



MODULE DBCSR_VAR_MOD

  USE dbcsr_types
  USE dbcsr_methods
  USE dbcsr_error_handling
  USE array_types,                     ONLY: array_data,&
       array_i1d_obj,&
       array_new,&
       array_nullify,&
       array_release,&
       array_size
  USE dbcsr_io
  USE dbcsr_operations
  USE dbcsr_ptr_util
  USE dbcsr_transformations
  USE dbcsr_util
  USE dbcsr_work_operations
  USE dbcsr_message_passing

  USE dbcsr_block_access
  USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
       dbcsr_iterator_next_block,&
       dbcsr_iterator_start,&
       dbcsr_iterator_stop

  USE dbcsr_dist_operations,           ONLY: create_bl_distribution,&
       dbcsr_get_stored_coordinates
  IMPLICIT NONE
  SAVE

  !*****************************************************************
  !sets up dbcsr/mpi variables

  TYPE(dbcsr_obj)                          :: matrix_a, matrix_b
  !labels as standard error, necessary dbcsr declaration
  TYPE(dbcsr_error_type)                   :: error

  !setting up variables for matrix and mpi
  INTEGER, DIMENSION(:), POINTER           :: rbs, cbs
  TYPE(array_i1d_obj)                      :: row_blk_sizes, col_blk_sizes, row_dist_a, col_dist_a
  INTEGER                                  :: npdims(2), ierr
  INTEGER, ALLOCATABLE, DIMENSION(:,:)     :: pgrid
  INTEGER prow, pcol, mp_comm, mp_group, nblkcols_total, nblkrows_total,proc_holds_blk
  TYPE(dbcsr_mp_obj)                       :: mp_env
  INTEGER                                  :: group, mynode, numnodes, myploc(2), N
  TYPE(dbcsr_distribution_obj)             :: dist_a, dist_b, dist_c
  REAL(LATTEPREC), DIMENSION(:), ALLOCATABLE       :: diag
  REAL(LATTEPREC)                                  :: my_block(2:2)
  LOGICAL                                  :: tr, found
  INTEGER, ALLOCATABLE, DIMENSION(:)       :: grid_dist
  TYPE(dbcsr_iterator)                     :: iter
  REAL(LATTEPREC)                                  :: chksum, chksum2
  INTEGER                                  :: TEMP

  !*****************************************************************
  !sets distribution to processors

CONTAINS

  SUBROUTINE myset_dist (dist_array, dist_size, nbins)
    TYPE(array_i1d_obj), INTENT(OUT)         :: dist_array
    INTEGER, INTENT(IN)                      :: dist_size, nbins

    INTEGER                                  :: i
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: grid_dist

    ALLOCATE (grid_dist(dist_size))
    CALL array_nullify (dist_array)

    FORALL (i = 1 : dist_size)
       grid_dist(i) = MODULO (nbins-i, nbins)
    END FORALL

    CALL array_new (dist_array, grid_dist, lb=1)
    DEALLOCATE (grid_dist)

  END SUBROUTINE myset_dist



  !*****************************************************************

END MODULE DBCSR_VAR_MOD

