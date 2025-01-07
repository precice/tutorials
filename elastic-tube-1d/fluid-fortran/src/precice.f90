module precice
  use, intrinsic :: iso_c_binding
  implicit none

  interface
  
    subroutine precicef_create(participantName, configFileName, &
      &                        solverProcessIndex, solverProcessSize, &
      &                        participantNameLength, configFileNameLength) &
      &  bind(c, name='precicef_create_')

      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: participantName
      character(kind=c_char), dimension(*) :: configFileName
      integer(kind=c_int) :: solverProcessIndex
      integer(kind=c_int) :: solverProcessSize
      integer(kind=c_int), value :: participantNameLength
      integer(kind=c_int), value :: configFileNameLength
    end subroutine precicef_create

    subroutine precicef_initialize() &
      &  bind(c, name='precicef_initialize_')

      use, intrinsic :: iso_c_binding
    end subroutine precicef_initialize

    subroutine precicef_advance(timestepLengthLimit) &
      &  bind(c, name='precicef_advance_')

      use, intrinsic :: iso_c_binding
      real(kind=c_double) :: timestepLengthLimit
    end subroutine precicef_advance

    subroutine precicef_finalize() &
      & bind(c, name='precicef_finalize_')

      use, intrinsic :: iso_c_binding
    end subroutine precicef_finalize

    subroutine precicef_requires_reading_checkpoint(isRequired) &
       &  bind(c, name='precicef_requires_reading_checkpoint_')

       use, intrinsic :: iso_c_binding
       integer(kind=c_int) :: isRequired
    end subroutine precicef_requires_reading_checkpoint

    subroutine precicef_requires_writing_checkpoint(isRequired) &
       &  bind(c, name='precicef_requires_writing_checkpoint_')

       use, intrinsic :: iso_c_binding
       integer(kind=c_int) :: isRequired
    end subroutine precicef_requires_writing_checkpoint

    subroutine precicef_get_mesh_dimensions(meshName, dimensions, &
      & meshNameLength) &
      & bind(c, name='precicef_get_mesh_dimensions_')

      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: meshName
      integer(kind=c_int) :: dimensions
      integer(kind=c_int), value :: meshNameLength
    end subroutine precicef_get_mesh_dimensions

    subroutine precicef_get_data_dimensions(meshName, dataName, &
      & dimensions, meshNameLength, dataNameLength) &
      & bind(c, name='precicef_get_data_dimensions_')

      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: meshName
      character(kind=c_char), dimension(*) :: dataName
      integer(kind=c_int) :: dimensions
      integer(kind=c_int), value :: meshNameLength
      integer(kind=c_int), value :: dataNameLength
    end subroutine precicef_get_data_dimensions

    subroutine precicef_is_coupling_ongoing(isOngoing) &
      &  bind(c, name='precicef_is_coupling_ongoing_')

      use, intrinsic :: iso_c_binding
      integer(kind=c_int) :: isOngoing
    end subroutine precicef_is_coupling_ongoing

    subroutine precicef_is_time_window_complete(isComplete) &
      &  bind(c, name='precicef_is_time_window_complete_')

      use, intrinsic :: iso_c_binding
      integer(kind=c_int) :: isComplete
    end subroutine precicef_is_time_window_complete     
    
    subroutine precicef_get_max_time_step_size(maxTimeStepSize) &
      & bind(c, name='precicef_get_max_time_step_size_')

      use, intrinsic :: iso_c_binding
      real(kind=c_double) :: maxTimeStepSize
    end subroutine precicef_get_max_time_step_size

    subroutine precicef_requires_mesh_connectivity_for(meshName, required, meshNameLength) &
      & bind(c, name='precicef_requires_mesh_connectivity_for_')

      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: meshName
      integer(kind=c_int) :: required
      integer(kind=c_int), value :: meshNameLength
    end subroutine precicef_requires_mesh_connectivity_for

    subroutine precicef_set_vertex(meshName, position, vertexID, meshNameLength) &
      &  bind(c, name='precicef_set_vertex_')

      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: meshName
      real(kind=c_double) :: coordinates(3)
      integer(kind=c_int) :: id
      integer(kind=c_int), value :: meshNameLength
    end subroutine precicef_set_vertex
    
    subroutine precicef_get_mesh_vertex_size(meshName, meshSize, meshNameLength) &
      &  bind(c, name='precicef_get_mesh_vertex_size_')

      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: meshName
      integer(kind=c_int) :: meshSize
      integer(kind=c_int), value :: meshNameLength
    end subroutine precicef_get_mesh_vertex_size    
    
    subroutine precicef_set_vertices(meshName, size, coordinates, ids, meshNameLength) &
      &  bind(c, name='precicef_set_vertices_')

      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: meshName
      integer(kind=c_int) :: size
      real(kind=c_double) :: coordinates(*)
      integer(kind=c_int) :: ids(*)
      integer(kind=c_int), value :: meshNameLength
    end subroutine precicef_set_vertices 

    subroutine precicef_set_edge(meshName, firstVertexID, secondVertexID, &
      &                          meshNameLength) &
      &  bind(c, name='precicef_set_edge_')

      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: meshName
      integer(kind=c_int) :: firstVertexID
      integer(kind=c_int) :: secondVertexID
      integer(kind=c_int), value :: meshNameLength
    end subroutine precicef_set_edge

    subroutine precicef_set_mesh_edges(meshName, size, ids, meshNameLength) &
      & bind(c, name='precicef_set_mesh_edges_')
    
      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: meshName
      integer(kind=c_int) :: size
      integer(kind=c_int) :: ids(*)
      integer(kind=c_int), value :: meshNameLength
    end subroutine precicef_set_mesh_edges

    subroutine precicef_set_triangle(meshName, firstEdgeID, secondEdgeID, &
      &                              thirdEdgeID, meshNameLength) &
      &  bind(c, name='precicef_set_triangle_')

      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: meshName
      integer(kind=c_int) :: firstEdgeID
      integer(kind=c_int) :: secondEdgeID
      integer(kind=c_int) :: thirdEdgeID
      integer(kind=c_int), value :: meshNameLength
    end subroutine precicef_set_triangle

    subroutine precicef_set_quad(meshName, firstVertexID, secondVertexID, &
      &                          thirdVertexID, fourthVertexID, &
      &                          meshNameLength ) &
      &  bind(c, name='precicef_set_quad_')

      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: meshName
      integer(kind=c_int) :: firstVertexID
      integer(kind=c_int) :: secondVertexID
      integer(kind=c_int) :: thirdVertexID
      integer(kind=c_int) :: fourthVertexID
      integer(kind=c_int), value :: meshNameLength
    end subroutine precicef_set_quad

    subroutine precicef_set_mesh_quads(meshName, size, ids, meshNameLength) &
      & bind(c, name='precicef_set_mesh_quads_')
    
      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: meshName
      integer(kind=c_int) :: size
      integer(kind=c_int) :: ids(*)
      integer(kind=c_int), value :: meshNameLength
    end subroutine precicef_set_mesh_quads

    subroutine precicef_set_tetrahedron(meshName, firstVertexID, secondVertexID, &
      &                                 thirdVertexID, fourthVertexID, &
      &                                 meshNameLength) &
      &  bind(c, name='precicef_set_tetrahedron_')

      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: meshName
      integer(kind=c_int) :: firstVertexID
      integer(kind=c_int) :: secondVertexID
      integer(kind=c_int) :: thirdVertexID
      integer(kind=c_int) :: fourthVertexID
      integer(kind=c_int), value :: meshNameLength
    end subroutine precicef_set_tetrahedron
  
    subroutine precicef_set_mesh_tetrahedra(meshName, size, ids, meshNameLength) &
      & bind(c, name='precicef_set_mesh_tetrahedra_')
    
      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: meshName
      integer(kind=c_int) :: size
      integer(kind=c_int) :: ids(*)
      integer(kind=c_int), value :: meshNameLength
    end subroutine precicef_set_mesh_tetrahedra

    subroutine precicef_requires_initial_data(isRequired) &
      & bind(c, name='precicef_requires_initial_data_')
    
      use, intrinsic :: iso_c_binding
      integer(kind=c_int) :: isRequired
    end subroutine precicef_requires_initial_data

    subroutine precicef_write_data(meshName, dataName, size, ids, &
      &                            values, meshNameLength, dataNameLength) &
      & bind(c, name='precicef_write_data_')
    
      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: meshName
      character(kind=c_char), dimension(*) :: dataName
      integer(kind=c_int) :: size
      integer(kind=c_int) :: ids(*)
      real(kind=c_double) :: values(*)
      integer(kind=c_int), value :: meshNameLength
      integer(kind=c_int), value :: dataNameLength
    end subroutine precicef_write_data

    subroutine precicef_read_data(meshName, dataName, size, ids, &
      &                           relativeReadTime, values, meshNameLength, &
      &                           dataNameLength) &
      & bind(c, name='precicef_read_data_')
    
      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: meshName
      character(kind=c_char), dimension(*) :: dataName
      integer(kind=c_int) :: size
      integer(kind=c_int) :: ids(*)
      real(kind=c_double) :: relativeReadTime
      real(kind=c_double) :: values(*)
      integer(kind=c_int), value :: meshNameLength
      integer(kind=c_int), value :: dataNameLength
    end subroutine precicef_read_data

    subroutine precicef_set_mesh_access_region(meshName, boundingBox, &
      &                                        meshNameLength) &
      & bind(c, name='precicef_set_mesh_access_region_')

      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: meshName
      real(kind=c_double) :: boundingBox(*)
      integer(kind=c_int), value :: meshNameLength
    end subroutine precicef_set_mesh_access_region

    subroutine precicef_get_mesh_vertex_ids_and_coordinates(meshName, size, &
      &                                                     ids, coordinates, &
      &                                                     meshNameLength) &
      & bind(c, name='precicef_get_mesh_vertex_ids_and_coordinates_')

      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: meshName
      integer(kind=c_int) :: size
      integer(kind=c_int) :: ids(*)
      real(kind=c_double) :: coordinates(*) 
      integer(kind=c_int), value :: meshNameLength
    end subroutine precicef_get_mesh_vertex_ids_and_coordinates

    subroutine precicef_requires_gradient_data_for(meshName, dataName, &
      &                                            required, meshNameLength, &
      &                                            dataNameLength) &
      & bind(c, name='precicef_requires_gradient_data_for_')

      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: meshName
      character(kind=c_char), dimension(*) :: dataName
      integer(kind=c_int) :: required
      integer(kind=c_int), value :: meshNameLength
      integer(kind=c_int), value :: dataNameLength
    end subroutine precicef_requires_gradient_data_for

    subroutine precicef_write_gradient_data(meshName, dataName, size, ids, &
      &                                     gradients, meshNameLength, &
      &                                     dataNameLength) &
      & bind(c, name='precicef_write_gradient_data_')

      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: meshName
      character(kind=c_char), dimension(*) :: dataName
      integer(kind=c_int) :: size
      integer(kind=c_int) :: ids(*)
      real(kind=c_double) :: gradients(*)
      integer(kind=c_int), value :: meshNameLength
      integer(kind=c_int), value :: dataNameLength
    end subroutine precicef_write_gradient_data

    subroutine precicef_get_version_information(versionInfo, lengthVersionInfo) &
      &        bind(c, name="precicef_get_version_information_")

      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: versionInfo
      integer(kind=c_int), value :: lengthVersionInfo
    end subroutine precicef_get_version_information

  end interface

end module precice
