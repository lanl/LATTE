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

SUBROUTINE SHUTDOWN_DBCSR

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


  CALL dbcsr_distribution_release(dist_a)
  CALL dbcsr_mp_release(mp_env)
  CALL array_release(row_dist_a)
  CALL array_release(col_dist_a)	
  CALL array_release(row_blk_sizes)
  CALL array_release(col_blk_sizes)


  CALL mp_world_finalize()


END SUBROUTINE SHUTDOWN_DBCSR
