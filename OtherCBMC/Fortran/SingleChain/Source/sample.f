      Subroutine Sample(Choice,R,W)
C      Samples the end-to-end distance of the chain
      Implicit None
      Include 'system.inc'

c     Global variables
      Integer Choice 
      Double Precision R,W
C    Choice:      parameter that specifies what the subroutine will do: 1 = initiallize; 2 =  sample; 3 = output final data 
C    R:           end-to-end distance of the current chain configuration
C    W:           weight of the current configuration.  The sampling scheme used in subroutine mcloop is called a dynamic scheme; in dynamics schemes W is always 1.  The subroutine Sample is written so it can also be used in static schemes where W is not always 1.

C     Local variables
      Integer I,Maxx
      Parameter(Maxx=501)     
      Double Precision normalize
      Double Precision histogram(Maxx)
      Double Precision dR
      Double Precision averageR2
C       Maxx:            number  of bins in the histogram    
C       normalize:       the variable name says it all...
C       histogram:       the variable name says it all...
C       dR:              inverse bin width in the histogram (converted to bin width when writing the output)
C       averageR2:       use it to calculate <R^2>   

C     Saving these variables between different calls to this subroutine
      Save normalize,histogram,dR,averageR2

C     If "Choice" = 1 then initialize variables for this subroutine
      If(Choice.Eq.1) Then
         Do I=1,Maxx
            histogram(I) = 0.0d0
         Enddo
         normalize = 0.0d0
         averageR2 = 0.0d0
         dR = Dble(Maxx-1)/Dble(Maxchain)
C        If the number of bins in the histogram is Maxx, why is dR not
C        Dble(Maxx)/Dble(Maxchain) ?
C        Answer: As we use function Int to populate the histogram, 
C        if dR=Dble(Maxx)/Dble(Maxchain) then a fully stretched chain of length
C        Maxchain  will fall outside the histogram and the program will crash.

C     If "Choice" = 2 then acumulate statistics on the chain end-to-end distances
      Elseif(Choice.Eq.2) Then
C        Int returns the largest integer whose absolute value does not exceed the absolute value of the argument and has the same sign as the argument 
         I      = 1 + Int(R*dR)
         histogram(I) = histogram(I) + W
         normalize    = normalize    + W

C        Start Modification to calculate R2: sum over all  R*R (with the correct weight). 
C           Use variable averageR2 to store the sum
         averageR2 =averageR2+R*R*W
C        End Modification to calculate R2

C     If "Choice" = 3 then output calculated quantities
      Elseif (Choice.Eq.3) Then
         normalize    = 1.0d0/normalize
         dR = 1.0d0/dR
         Write(23,*) "# R (sigma)              P(R)/dR" 

         Do I=1,Maxx
            if (histogram(I) .GT. 0.5) THEN 
              Write(23,*) Dble(I-1)*dR+dR/2.,histogram(I)*
     &normalize/dR
            ENDIF
         Enddo
C        Start Modification to calculate R2
C          Normalize and print R2
      print*,"averageR2         : ",averageR2*normalize
C        End Modification to calculate R2
         
      Endif
      Return
      End
