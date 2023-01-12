      Subroutine Ran_Sphere(X,Y,Z)
C       Random number from a uniform distribution on a sphere: 
C       returns the X, Y and Z coordinates of a point on a sphere of radius 1 and center (0,0,0)
C       This subroutine follows the algorithm proposed by Marsaglia (1972):
C       Marsaglia, G. "Choosing a Point from the Surface of a Sphere." Ann. Math. Stat. 43, 645-646, 1972. 
      Implicit None

C     Global variables
      Double Precision X,Y,Z

C     Local variables
      Double Precision Ran1,Ran2,Ranh,Ransq,RandomNumber

C      Note: RandomNumber() is a call to a function that returns a random number from a uniform distribution in the interval [0,1] 
 100  Ran1  = 2.0d0*RandomNumber() - 1.0d0
      Ran2  = 2.0d0*RandomNumber() - 1.0d0
      Ransq = Ran1*Ran1 + Ran2*Ran2
 
      If (Ransq.Ge.1.0d0) Goto 100
 
      Ranh = 2.0d0*Dsqrt(1.0d0 - Ransq)
      X = Ran1*Ranh
      Y = Ran2*Ranh
      Z = 1.0d0 - 2.0d0*Ransq
 
      Return
      End
