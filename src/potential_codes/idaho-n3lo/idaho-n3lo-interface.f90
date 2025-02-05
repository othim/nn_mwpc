module idaho_n3lo_interface

  implicit none
  
contains

  subroutine vidaho_n3lo(pout, pin, coup_, S, j_, T, tz, POT) &
          bind(C,name="idaho_n3lo_fort_interface")

    REAL(8) ,INTENT(IN)    :: pout 
    REAL(8) ,INTENT(IN)    :: pin  
    LOGICAL,INTENT(IN)    :: coup_  
    INTEGER,INTENT(IN)    :: S      
    INTEGER,INTENT(IN)    :: j_     
    INTEGER,INTENT(IN)    :: T   ! not used/tested in interface. only in pnijm   
    INTEGER,INTENT(IN)    :: tz     
    REAL(8),INTENT(INOUT) :: POT(6) 
    
    INTEGER  :: inn, j, kread, kwrite, kda, kpunch                                                   
    REAL(8)  :: v,xmev,ymev                                                     
    CHARACTER (LEN=4) :: label                                                  
    LOGICAL :: coup, sing, trip, heform, endep     
    
    ! Interaction commion block
    COMMON /crdwrt/ kread,kwrite,kpunch,kda(9)
    ! arguments and values of this subroutine:
    COMMON /cpot/   v(6),xmev,ymev
    COMMON /cstate/ j,heform,sing,trip,coup,endep,label
    COMMON /cnn/ inn

    
    coup = .FALSE. ; sing = .FALSE. ; trip = .FALSE. ; heform = .FALSE.
    
    SELECT CASE(tz)
    CASE(-1)
       ! pp
       inn = 1
    CASE(0)
       ! pn
       inn = 2
    CASE(+1)
       ! nn
       inn = 3
    END SELECT

    IF (coup_)  coup = .TRUE.
    IF (S == 1) trip = .TRUE.
    IF (S == 0) sing = .TRUE.

    j = j_

    kread  = 5
    kwrite = 6 
    
    xmev = pout
    ymev = pin
    
    call n3lo
    POT = v

    end subroutine vidaho_n3lo

end module idaho_n3lo_interface
