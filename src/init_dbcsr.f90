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

SUBROUTINE INIT_DBCSR

  USE DBCSR_VAR_MOD
  USE dbcsr_config
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

  USE CONSTANTS_MOD

  IMPLICIT NONE

  !sets mpi

  !sets up dbcsr matrix
  ! the matrix will contain nblkrows_total row blocks and nblkcols_total column blocks


  !initiallizing mpi

  CALL mp_world_init(mp_comm)

  npdims(:) = 0

  CALL mp_cart_create (mp_comm, 2, npdims, myploc, group)	

  CALL mp_environ (numnodes, mynode, group)

  ALLOCATE (pgrid(0:npdims(1)-1, 0:npdims(2)-1))

  DO prow = 0, npdims(1)-1
     DO pcol = 0, npdims(2)-1
        CALL mp_cart_rank (group, (/ prow, pcol /), pgrid(prow, pcol))
     ENDDO
  ENDDO

  ! Create the dbcsr_mp_obj
  CALL dbcsr_mp_new (mp_env, pgrid, group, mynode, numnodes,&
       myprow=myploc(1), mypcol=myploc(2))

  DEALLOCATE(pgrid)

  ! Use BLAS rather than the SMM

  CALL dbcsr_set_conf_mm_driver(2, error=error)

  ! Now with padding

  nblkrows_total=(HDIM-1)/BLKSZ + BLKSZ
  nblkcols_total=(HDIM-1)/BLKSZ + BLKSZ


  !sets the block size for each row and column
  ALLOCATE(rbs(nblkrows_total))
  ALLOCATE(cbs(nblkcols_total))  
  rbs(:)=BLKSZ
  cbs(:)=BLKSZ

  CALL array_nullify (row_blk_sizes)
  CALL array_nullify (col_blk_sizes)
  CALL array_new (row_blk_sizes, rbs, gift=.TRUE.)
  CALL array_new (col_blk_sizes, cbs, gift=.TRUE.)


  !sets distribution to processors
  CALL myset_dist (row_dist_a, nblkrows_total, npdims(1))
  CALL myset_dist (col_dist_a, nblkcols_total, npdims(2))


  !Sets the distribution object
  CALL dbcsr_distribution_new (dist_a, mp_env, row_dist_a, col_dist_a)

END SUBROUTINE INIT_DBCSR
