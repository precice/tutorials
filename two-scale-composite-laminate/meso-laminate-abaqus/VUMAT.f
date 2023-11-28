      include './get_strains.f'
      include './material_model.f'
! Authors: Ishaan Desai (desaii)
      subroutine vumat(
! Read only (unmodifiable)variables -
     1     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2     stepTime, totalTime, dt, cmname, coordMp, charLength,
     3     props, density, strainInc, relSpinInc,
     4     tempOld, stretchOld, defgradOld, fieldOld,
     5     stressOld, stateOld, enerInternOld, enerInelasOld,
     6     tempNew, stretchNew, defgradNew, fieldNew,
! Write only (modifiable) variables -
     7     stressNew, stateNew, enerInternNew, enerInelasNew)
!
      include 'vaba_param.inc'

      dimension props(nprops), density(nblock), coordMp(nblock, *),
     1     charLength(nblock), strainInc(nblock, ndir + nshr),
     2     relSpinInc(nblock, nshr), tempOld(nblock),
     3     stretchOld(nblock, ndir + nshr),
     4     defgradOld(nblock, ndir + nshr + nshr),
     5     fieldOld(nblock, nfieldv), stressOld(nblock, ndir + nshr),
     6     stateOld(nblock, nstatev), enerInternOld(nblock),
     7     enerInelasOld(nblock), tempNew(nblock),
     8     stretchNew(nblock, ndir + nshr),
     8     defgradNew(nblock, ndir + nshr + nshr),
     9     fieldNew(nblock, nfieldv),
     1     stressNew(nblock, ndir + nshr), stateNew(nblock, nstatev),
     2     enerInternNew(nblock), enerInelasNew(nblock)
!
      character*80 cmname

      real*8, dimension(nblock, nstatev) :: state
      real*8, dimension(nblock, ndir + nshr) :: strains, stresses
      integer :: d, counter

      ! preCICE variables
      integer :: rank, size, ongoing, dimensions, bool
      double precision :: preCICE_dt, abaqus_dt ! We cannot change the dt variable
      double precision, dimension(nblock*ndir) :: strains1to3, strains4to6,
     * stresses1to3, stresses4to6, couplingVertices
      integer, dimension(nblock) :: vertexIDs

      write(*,*) "t = ", totalTime, ",", " dt = ", dt, " VUMAT: Entering VUMAT."

      ! When not in initialization, vertexIDs of preCICE are defined
      ! manually, because VUMAT cannot hold global data.
      ! TODO: Need to find a way to store global data in VUMAT
      do k = 0, nblock - 1
         vertexIDs(k + 1) = k
      end do

      if (dt /= 1) then
         if (totalTime < 2.*dt) then ! only the first step

            ! Get MPI rank and size (total number of MPI processors in this job)
            call vgetnumcpus(size)
            call vgetrank(rank)

            ! Create preCICE participant
            call precicef_create("Laminate-3D-ply",
     *       "/home/desaii/composite-multiscale/precice-config.xml",
     *       rank, size)

            ! Get problem dimensions from preCICE
            call precicef_get_mesh_dimensions("laminate-meso-mesh", dimensions)

            counter = 1
            do k = 1, nblock
               do d = 1, ndir
                  couplingVertices(counter) = coordMp(k, d)
                  counter = counter + 1
               end do
            end do

            ! Set coupling mesh vertices in preCICE
            call precicef_set_vertices("laminate-meso-mesh",
     *       nblock, couplingVertices, vertexIDs)

            call get_strains(nblock, ndir, nshr, nstatev,
     *       strainInc, stateOld, state, strains)

            ! Apply the calculated stresses
            do k = 1, nblock
               stateNew(k, :) = state(k, :)
            end do

            counter = 1
            do k = 1, nblock

               do d = 1, ndir
                  strains1to3(counter) = strains(k, d)
                  strains4to6(counter) = strains(k, d+ndir)
                  counter = counter + 1

               end do ! ndir

               do d = 1, nstatev

                  stateNew(k, d) = state(k, d)
               
               end do ! nstatev

            end do ! nblock

            call precicef_requires_initial_data(bool)

            if (bool == 1) then
               call precicef_write_data("laminate-meso-mesh",
     *          "strains1to3", nblock, vertexIDs, strains1to3)
               call precicef_write_data("laminate-meso-mesh",
     *          "strains4to6", nblock, vertexIDs, strains4to6)
            end if

            write(*,*) "VUMAT: Initial data written to preCICE."

            ! Initialize preCICE
            call precicef_initialize()
            call precicef_get_max_time_step_size(preCICE_dt)

            write(*,*) "VUMAT: preCICE initialization complete."

         end if ! if (totalTime < 2.*dt)

         ! Check if coupling is still going on
         call precicef_is_coupling_ongoing(bool)
         !call assert(bool.eq.1)

! Get strains and write them to preCICE ===============================
         call get_strains(nblock, ndir, nshr, nstatev,
     *    strainInc, stateOld, state, strains)
 
         counter = 1
         do k = 1, nblock

            do d = 1, ndir

               strains1to3(counter) = strains(k, d)
               strains4to6(counter) = strains(k, d+ndir)
               counter = counter + 1

            end do ! ndir

            do d = 1, nstatev
               stateNew(k, d) = state(k, d)
            end do ! nstatev

         end do ! nblock

         call precicef_write_data("laminate-meso-mesh",
     *       "strains1to3", nblock, vertexIDs, strains1to3)
         call precicef_write_data("laminate-meso-mesh",
     *       "strains4to6", nblock, vertexIDs, strains4to6)

         write(*,*) "(t = ", totalTime, ") VUMAT: Strains written to preCICE."

! ==========================================================================

         call precicef_get_max_time_step_size(preCICE_dt)
         abaqus_dt = min(preCICE_dt, dt)

         call precicef_advance(abaqus_dt)

         write(*,*) "(t = ", totalTime, ") VUMAT: Coupling has been advanced."

! Read stresses from preCICE and apply them ===============================

         if (totalTime > dt) then ! from the second step onward

            call precicef_read_data("laminate-meso-mesh",
     *       "stresses1to3", nblock, vertexIDs, dt, stresses1to3)
            call precicef_read_data("laminate-meso-mesh",
     *       "stresses4to6", nblock, vertexIDs, dt, stresses4to6)

            ! Loop through material points to apply stresses
            counter = 1
            do k = 1, nblock
               do d = 1, ndir
                  stressNew(k, d) = stresses1to3(counter)
                  stressNew(k, d+ndir) = stresses4to6(counter)
                  counter = counter + 1

               end do ! ndir
            end do ! nblock

            !write(*,*) "stressNew: ", stressNew

            write(*,*) "(t = ", totalTime, ") VUMAT: Stresses applied."

         else ! only in the first step
            call solve_material_model(nblock, ndir, nshr,
     *       nprops, props, strains, stresses) 

            ! Apply the calculated stresses
            do k = 1, nblock
               stressNew(k, :) = stresses(k, :)
            end do

            write(*,*) "(t = ", totalTime, ") VUMAT: stresses not read from preCICE, but instead the material model solved."

         end if
! ==========================================================================

         if (totalTime == 1.0) then ! Hardcode the end time because we cannot access it
            call precicef_finalize()
         end if

         write(*,*) "(t = ", totalTime, ") VUMAT: run complete."

      else ! if (dt == 1)

         call get_strains(nblock, ndir, nshr, nstatev,
     *    strainInc, stateOld, state, strains)

         ! Apply the calculated stresses
         do k = 1, nblock
            stateNew(k, :) = state(k, :)
         end do

         call solve_material_model(nblock, ndir, nshr,
     *    nprops, props, strains, stresses)

         ! Apply the calculated stresses
         do k = 1, nblock
            stressNew(k, :) = stresses(k, :)
         end do

         write(*,*) "(t = ", totalTime, ") VUMAT: Material model solved at dt = 1"

      end if ! if (dt /= 1)

      return
      end ! Subroutine vumat
