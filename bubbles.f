!----------------------------------------------------------------------
!       Module for shared bubble variables or arrays
!----------------------------------------------------------------------
      module bubbleGen

        integer nbub
        real*8  rbub
        real*8, allocatable:: bubcoord(:,:)

      end module
!----------------------------------------------------------------------
        PROGRAM bubbles
        use bubbleGen
        IMPLICIT NONE
!----------------------------------------------------------------------
!       This code is specialized for single subchannel geometry with
!       spacer grid and mixing vanes.
!
! This routine is developed to created the required bubble information
! such as, bubble ID, bubble coordiantes and radii. And the Monta Carlo
! rejection sampling method is used to avoid the intersection between
! bubbles and the walls or between different bubbles
!
! Some variables:
! mixvanes(2)           The range of mixing vanes region
! coordrep              The flag for coordinate replacement
!
! Jun Fang, 2015.
!----------------------------------------------------------------------

        INTEGER  ierror,        i,j,k
        integer::clock,   n_dimension=4, size
        integer  ibub
        REAL*8   rbox1(6),  rbox2(3), mixvanes(2)
        real*8   distance, mind2w, mindist,  mindistmp, r
        real*8   ghostupplim, ghostlowlim
        real*8   random_replace(4),  temp(1,3), distw(4)
        REAL*8, ALLOCATABLE:: minw(:),
     &                        random1d(:),   random2d(:,:)
        LOGICAL::gate=.TRUE., check,         gothroughlist,
     &           coordrep=.FALSE.

        INTEGER::nghost,ntotal
        REAL*8:: ghost_area_ratio,outline(3),bubtmp(3)
        REAL*8, ALLOCATABLE:: bublist(:,:)
        real*8, allocatable:: bubdupl(:,:)        

        CHARACTER (LEN=3) :: DomainType
!----------------------------------------------------------------------
!       Creating the output files 
!       Reading input parameters from files (options.dat)
!----------------------------------------------------------------------
        
        OPEN(UNIT=1,FILE='options.dat',IOSTAT=ierror)
        IF(ierror .NE. 0)STOP"Error opening file options.txt"

        OPEN(UNIT=2,FILE='bubbles.inp',IOSTAT=ierror)
        IF(ierror .NE. 0)STOP"Error opening file output.txt"

        OPEN(UNIT=6,FILE='Logs.out',IOSTAT=ierror)
        IF(ierror .NE. 0)STOP"Error opening file output.txt"

!        OPEN(unit=5,FILE='void_x.dat',IOSTAT=ierror)
!        IF(ierror .NE. 0)STOP"Error opening file void_x.dat"

        READ(1,*)
        READ(1,*)DomainType             !Domain type
        READ(1,*)
        READ(1,*)nbub                   !number of bubbles
        READ(1,*)
        READ(1,*)rbub                   !bubble radius
        READ(1,*)
        READ(1,*)rbox1(:)               !domain dimensions
        READ(1,*)
        READ(1,*)mixvanes(:)            !range of mixing vanes
        READ(1,*)
        READ(1,*)r                      !rod radius
        READ(1,*)
        READ(1,*)mindist                !min distance b/w bubbles (in the unit of bubble radii)
        READ(1,*)
        READ(1,*)ghost_area_ratio       !ghost region ratio

        mindist = mindist * rbub        !min distance 

        !write(*,*) 'Min. interbubble distance =', mindist
        if (DomainType .eq. "Sub") then
           write(*,*) 'Domain Type: PWR subchannel'
           write(*,*) 'Number of bubbles =', nbub
           write(*,*) 'Bubble radius =', rbub
           if(mixvanes(1).eq.mixvanes(2)) then
           write(*,*) 'Consider Suchannel alone!'
           else
           write(*,*) 'Consider SGMV region!'
           write(*,"(A, ES16.6, A, ES16.6)")
     &                ' SGMV region ranges from', mixvanes(1),
     &                ' to', mixvanes(2)
           endif
        endif

!----------------------------------------------------------------------
!       The initialization of bubble center coordinates
!----------------------------------------------------------------------
!Subtracting the radius from the box walls, gives room, lets me not worry
! about crossing at wall
        DO i=2,6,2
           rbox2(i/2)=rbox1(i)-rbox1(i-1)
        END DO
        outline(:)=rbox2(:)

        ghost_area_ratio=ghost_area_ratio/100.0

!Forming the random number arrays and the bubble coordinate arrays.
        allocate( random1d(nbub*3), random2d(nbub,3) )  
        ALLOCATE( bubcoord(nbub,3), minw(nbub) )

!forming the random array
        CALL SYSTEM_CLOCK(COUNT=clock)
        CALL random_number_init(nbub*3,clock)
        CALL RANDOM_NUMBER(random1d)

        DO i=1,nbub
           random2d(i,:)=random1d((i*3-2):i*3)
        END DO

!       Initialize the bubble center coordinates
        DO i=1,nbub
           DO j=1,3
              bubcoord(i,j)=random2d(i,j)*rbox2(j) !+ rbox1(2*j-1)
           END DO
        END DO

!       Determine the bubble positions in SGMV region first, and outside
!       the SGMV, bubbles will be initialized randomly. 
        ibub = 0        !The number of determined bubbles
        if(mixvanes(1).lt.mixvanes(2)) then     !put two identical numbers to exclude auxiliary region
           call uniformBubble(ibub, r, rbox1, mixvanes)
        endif
        write(*,*) 'number of uniform bubbles =', ibub
!----------------------------------------------------------------------
!       check if the first bubble center is in the right place, replace it if
!       not
!       The radius of fuel rod (r) is increased a little to avoid that
!       bubble's initial position is too close to walls
!----------------------------------------------------------------------
        i = ibub + 1
	call d2wall(mind2w, bubcoord(i,1:3), rbox2, r, mixvanes)

!...Check if the bubble is in SGMV region
        call checkregion(bubcoord(i,1), coordrep, mixvanes, rbub)
        IF(mind2w.LE.(2.0d0*rbub) .or. coordrep) THEN
          L0:DO
          IF(gate)THEN
           CALL random_number_init(n_dimension,int(random1d(i)*1000000))
           CALL RANDOM_NUMBER(random_replace)
           gate=.FALSE.
          ELSE
           CALL random_number_init( n_dimension, 
     &                          int(random_replace(2)*1000000))
           CALL RANDOM_NUMBER(random_replace)
          END IF

          DO k=1,3
           bubcoord(i,k)=rbox2(k)*random_replace(k+1) !+ rbox1(2*k-1)
          END DO
!...Check if the bubble is in SGMV region
          call checkregion(bubcoord(i,1), coordrep, mixvanes, rbub)
!...  Check the wall again !
          call d2wall(mind2w, bubcoord(i,1:3), rbox2, r, mixvanes)
          IF(mind2w.le.(2.0d0*rbub) .or. coordrep) CYCLE L0

          EXIT L0
          END DO L0
        ENDIF
        write(*,*) 'The first random bubble center coordinates: '
        write(*,"(3(1xES16.6))") bubcoord(i,:)
        write(*,*)

!----------------------------------------------------------------------
!       Determine the rest bubble centers one by one
!----------------------------------------------------------------------
        DO i=ibub+2,nbub
!...Check if the bubble is in mixing vanes region
           call checkregion(bubcoord(i,1), coordrep, mixvanes, rbub)

!...Find the minimum distance to the walls
           call d2wall(mind2w, bubcoord(i,1:3), rbox2, r, mixvanes)

!...Find the minimum distance to other bubbles
           mindistmp = 1e9
           do j=ibub+1,i-1  !loop over the centers that have already been found
              distance=sqrt(((bubcoord(i,1)-bubcoord(j,1))**2)+
     &                      ((bubcoord(i,2)-bubcoord(j,2))**2)+
     &                      ((bubcoord(i,3)-bubcoord(j,3))**2))
              if(mindistmp.gt.distance) mindistmp = distance
           end do

           gate = .TRUE.
!...Check the overlap
        if ((mindistmp .LE. mindist) .OR. (mind2w.LE.(2.0d0*rbub)) 
     &       .or. coordrep) then
           L1: do
           IF(gate)THEN
              CALL random_number_init(n_dimension,
     &                          int(random1d(i)*1000000))
             CALL RANDOM_NUMBER(random_replace)
              gate=.FALSE.
           ELSE
              CALL random_number_init(n_dimension,
     &                          int(random_replace(2)*1000000))
              CALL RANDOM_NUMBER(random_replace)
           END IF
           DO k=1,3
              bubcoord(i,k)=rbox2(k)*random_replace(k+1) !+ rbox1(2*k-1)
           END DO
!...Check if the bubble is in mixing vanes region
           call checkregion(bubcoord(i,1), coordrep, mixvanes, rbub)
!...Find the minimum distance to the walls
           call d2wall(mind2w, bubcoord(i,1:3), rbox2, r, mixvanes)
!...Find the minimum distance to other bubbles
           mindistmp = 1e9
           do j=ibub+1,i-1  !loop over the centers that have already been found
              distance=sqrt(((bubcoord(i,1)-bubcoord(j,1))**2)+
     &                      ((bubcoord(i,2)-bubcoord(j,2))**2)+
     &                      ((bubcoord(i,3)-bubcoord(j,3))**2))
              if(mindistmp.gt.distance) mindistmp = distance
           end do
           if ((mindistmp .LE. mindist).OR.(mind2w.LE.(2.0d0*rbub))
     &         .or. coordrep) cycle L1
           exit L1
           end do L1
        end if

        END DO  !for i loop

!----------------------------------------------------------------------
!       Post-processing
!----------------------------------------------------------------------
        DO i=ibub+1,nbub   
           DO j=1,3
              IF(j .EQ. 1)THEN
                 bubcoord(i,j)=bubcoord(i,j)+rbox1(1)
              ELSE IF(j .EQ. 2)THEN
                 bubcoord(i,j)=bubcoord(i,j)+rbox1(3)
              ELSE IF(j .EQ. 3)THEN
                 bubcoord(i,j)=bubcoord(i,j)+rbox1(5)
              END IF
           END DO
        END DO

!...Double check the coordinates we got 
        check=.false.
        DO i=1,nbub
           distw(1)=sqrt((bubcoord(i,3)-rbox1(5))**2 + 
     &                  (bubcoord(i,2)-rbox1(3))**2)
           distw(2)=sqrt((bubcoord(i,3)-rbox1(6))**2 + 
     &                  (bubcoord(i,2)-rbox1(3))**2)
           distw(3)=sqrt((bubcoord(i,3)-rbox1(5))**2 + 
     &                  (bubcoord(i,2)-rbox1(4))**2)
           distw(4)=sqrt((bubcoord(i,3)-rbox1(6))**2 + 
     &                  (bubcoord(i,2)-rbox1(4))**2)
           minw(i)=distw(1)
           DO k=2,4
              IF(minw(i).gt.distw(k)) minw(i)=distw(k)
           ENDDO

        IF(minw(i).lt.(r+2.0d0*rbub))THEN
           check=.true.
           WRITE(*,*)
           WRITE(*,*) 'Bubble touches the wall!!!'
           WRITE(*,*) 'The problem point is',bubcoord(i,:)
           WRITE(*,*) 'Distance is',(minw(i)-(r))
           WRITE(*,*)
        END IF

        DO j=1,nbub
           IF(i.EQ.j)CYCLE
           distance=sqrt(((bubcoord(i,1)-bubcoord(j,1))**2)+
     &                   ((bubcoord(i,2)-bubcoord(j,2))**2)+
     &                   ((bubcoord(i,3)-bubcoord(j,3))**2))
           IF(distance.lt.mindist)THEN
              check=.true.
              WRITE(*,*) 'Bubble touches another bubble!!!'
           END IF
        END DO

        ENDDO

        if(check) then
           WRITE(*,*)"Problem was found!!!"
        else
           WRITE(*,*)"So far, so good!!!"
        endif
!----------------------------------------------------------------------
!       Bubble sorting (based on x coordinates) 
!----------------------------------------------------------------------
        size=nbub
        gothroughlist=.TRUE.
        DO WHILE (gothroughlist)
                gothroughlist=.FALSE.
                DO i=1,size-1
                        IF (bubcoord(i,1) .GT. bubcoord(i+1,1)) THEN
                                temp(1,:)=bubcoord(i,:)
                                bubcoord(i,:)=bubcoord(i+1,:)
                                bubcoord(i+1,:)=temp(1,:)
                                gothroughlist=.TRUE.
                        END IF
                END DO
                size=size-1
        END DO

!...Find the ghost bubbles 
        WRITE(*,"(A, 6(1xES16.6))") ' Dimensions: ',rbox1(:)
        WRITE(*,"(A, f5.2, A)") ' The size of ghost domain is ',
     &                  ghost_area_ratio*100,'%'
        nghost=0
        DO i=1,nbub
           DO j=1,3
              ghostupplim = rbox1(2*j)  -ghost_area_ratio*outline(j)
              ghostlowlim = rbox1(2*j-1)+ghost_area_ratio*outline(j)
              IF (bubcoord(i,j) .GT. ghostupplim) then
                  nghost=nghost+1
              ELSEIF (bubcoord(i,j) .LT. ghostlowlim) then
                  nghost=nghost+1
              ENDIF
           END DO
        END DO

        WRITE(*,*) nghost, 'ghost bubbles generated!!!'
        ntotal=nbub+nghost
        ALLOCATE (bublist(ntotal,4))
        bublist(:,:) = 0.0

        DO i=1,nbub
           bublist(i,1) = REAL(i)
           DO j=1,3
              bublist(i,j+1)=bubcoord(i,j)
           ENDDO
        ENDDO

        nghost = 1
        DO i=1,nbub
           DO j=1,3
              ghostupplim = rbox1(2*j)  -ghost_area_ratio*outline(j)
              ghostlowlim = rbox1(2*j-1)+ghost_area_ratio*outline(j)
              IF (bubcoord(i,j) .GT. ghostupplim) THEN
                 bublist(nbub+nghost,2)=bubcoord(i,1)
                 bublist(nbub+nghost,3)=bubcoord(i,2)
                 bublist(nbub+nghost,4)=bubcoord(i,3)
                 bublist(nbub+nghost,j+1)=bubcoord(i,j)-outline(j)
                 WRITE(6,*) 'Bubble',i,'has the ghost bubble'
                 nghost=nghost+1
              ELSEIF (bubcoord(i,j) .LT. ghostlowlim) THEN
                 bublist(nbub+nghost,2)=bubcoord(i,1)
                 bublist(nbub+nghost,3)=bubcoord(i,2)
                 bublist(nbub+nghost,4)=bubcoord(i,3)
                 bublist(nbub+nghost,j+1)=bubcoord(i,j)+outline(j)
                 WRITE(6,*) 'Bubble',i,'has the ghost bubble'
                 nghost=nghost+1
              ENDIF
           ENDDO
        ENDDO

        DO i=1,ntotal
           WRITE(2,1010) int(bublist(i,1)), bublist(i,2:4), rbub
1010    FORMAT(I8, ES16.6, ES16.6, ES16.6, ES16.6)
        END DO

        deallocate( random1d, random2d )
        deallocate( bubcoord, minw )


      END PROGRAM


!----------------------------------------------------------------------
!       The subroutine for random seed 
!----------------------------------------------------------------------

      SUBROUTINE random_number_init(n,input)
        IMPLICIT NONE
        INTEGER::i,n,input
        INTEGER, ALLOCATABLE::seed(:)

        CALL RANDOM_SEED(size=n)
        ALLOCATE(seed(n))

        seed=input+37*(/(i-1,i=1,n)/)
        CALL RANDOM_SEED(PUT=seed)

        DEALLOCATE(seed)
      END SUBROUTINE


!----------------------------------------------------------------------
!       The subroutine to check if a bubble is in mixing vanes or not
!----------------------------------------------------------------------
      subroutine checkregion(xpostn, coordrep, mixvanes, rbub)
        implicit none
        logical :: coordrep
        real*8     xpostn, mixvanes(2)
        real*8     rbub

        if(xpostn.ge.(mixvanes(1)-rbub) .and.
     &     xpostn.le.(mixvanes(2)+rbub) ) then
           coordrep = .True.
        elseif(xpostn.ge.0.0405-rbub) then
!       The hard-coded parameter prevents initial bubbles across outlet
           coordrep = .True.
        else
!           write(*,*) xpostn, mixvanes
           coordrep = .False.
        endif
        return
      end subroutine


!----------------------------------------------------------------------
!       The subroutine to find out the minimum bubble wall distance
!----------------------------------------------------------------------
      subroutine d2wall(mind2w, bubcoordtmp, rbox2, r, mixvanes)
	implicit none
        integer j
        real*8 r
	real*8 bubcoordtmp(1,3), distw(12), mind2w, rbox2(3)
        real*8 mixvanes(2)
		
!       Four rod surfaces
        distw(1)=sqrt((bubcoordtmp(1,3))**2 + (bubcoordtmp(1,2))**2)
        distw(2)=sqrt((bubcoordtmp(1,3)-rbox2(3))**2 +
     &                (bubcoordtmp(1,2))**2)
        distw(3)=sqrt((bubcoordtmp(1,3))**2 + 
     &                (bubcoordtmp(1,2)-rbox2(2))**2)
        distw(4)=sqrt((bubcoordtmp(1,3)-rbox2(3))**2 + 
     &                (bubcoordtmp(1,2)-rbox2(2))**2)
        distw(1:4) = distw(1:4) - r 
!       Inlet/Outlet faces
        distw(5)=abs(bubcoordtmp(1,1))
        distw(6)=abs(bubcoordtmp(1,1)-rbox2(1))
!       Periodic faces
        distw(7)=abs(bubcoordtmp(1,2))
        distw(8)=abs(bubcoordtmp(1,2)-rbox2(2))        
        distw(9)=abs(bubcoordtmp(1,3))
        distw(10)=abs(bubcoordtmp(1,3)-rbox2(3))

!       The beginning and ending of spacer grid and mixing vanes region
        distw(11)=abs(bubcoordtmp(1,1)-mixvanes(1))
        distw(12)=abs(bubcoordtmp(1,1)-mixvanes(2))

        mind2w  =distw(1)
        DO j = 2, 12
           IF(mind2w.GE. distw(j)) mind2w=distw(j)
        END DO
	return
      end subroutine

!----------------------------------------------------------------------
!       The subroutine to generate bubble information in the spacer grid and
!       mixing vanes region
!----------------------------------------------------------------------
      subroutine uniformBubble(ibub, r, rbox1, mixvanes)
        use bubbleGen
        implicit none

        integer  ibub, i
        real*8   x1, x2, x3, cl(2), d2w, d2w_center
        real*8   pi, d1, d2, theta1, theta2, theta3
        real*8   rbox1(6), mixvanes(2)
        real*8   r, delta

        cl(1)   = (rbox1(3)+rbox1(4))/2.0d0
        cl(2)   = (rbox1(5)+rbox1(6))/2.0d0
        d2w_center    = sqrt((rbox1(4)-rbox1(3))**2.0d0 +
     &                       (rbox1(6)-rbox1(5))**2.0d0)/2.0d0 - r

        pi     = 4.0d0*atan(1.0d0)
        x1     = mixvanes(1) + 3.0d0*rbub
        d1     = d2w_center/2.0
        d2     = d1 + 3.0*rbub
        delta  = 5.1d0                  !in the unit of bubble radius
        theta1 = pi*(45.0/180.0)
        theta2 = pi*(25.0/180.0)
        theta3 = pi*(65.0/180.0)

        do while (x1.lt.0.012)        !the x limit is hardcoded for
                                      !reduced size subchannel
           do i = 1,4
              bubcoord(ibub+i,1) = x1
              bubcoord(ibub+i,2) = cl(1) +
     &               d1*cos(theta1+real(i-1)*pi/2.0d0)
              bubcoord(ibub+i,3) = cl(2) +
     &               d1*sin(theta1+real(i-1)*pi/2.0d0)
              bubcoord(ibub+4+i,1) = x1
              bubcoord(ibub+4+i,2) = cl(1) +
     &               d2*cos(theta2+real(i-1)*pi/2.0d0)
              bubcoord(ibub+4+i,3) = cl(2) +
     &               d2*sin(theta2+real(i-1)*pi/2.0d0)
              bubcoord(ibub+8+i,1) = x1
              bubcoord(ibub+8+i,2) = cl(1) +
     &               d2*cos(theta3+real(i-1)*pi/2.0d0)
              bubcoord(ibub+8+i,3) = cl(2) +
     &               d2*sin(theta3+real(i-1)*pi/2.0d0)
           enddo
           x1   = x1 + delta*rbub     !This multiplier is hardcoded
           ibub = ibub + 12
        enddo
        do while(x1.gt.0.012.and.x1.lt.0.01717)
           do i = 1,2
              bubcoord(ibub+i,1) = x1
              bubcoord(ibub+i,2) = cl(1) +
     &               d1*cos(theta1+real(i-1)*pi)
              bubcoord(ibub+i,3) = cl(2) +
     &               d1*sin(theta1+real(i-1)*pi)
              bubcoord(ibub+2+i,1) = x1
              bubcoord(ibub+2+i,2) = cl(1) +
     &               d2*cos(theta2+real(i-1)*pi)
              bubcoord(ibub+2+i,3) = cl(2) +
     &               d2*sin(theta2+real(i-1)*pi)
              bubcoord(ibub+4+i,1) = x1
              bubcoord(ibub+4+i,2) = cl(1) +
     &               d2*cos(theta3+real(i-1)*pi)
              bubcoord(ibub+4+i,3) = cl(2) +
     &               d2*sin(theta3+real(i-1)*pi)
           enddo
           x1   = x1 + delta*rbub     !This multiplier is hardcoded
           ibub = ibub + 6
        enddo

        return
      end subroutine
