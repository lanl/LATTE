subroutine LoadPairPotParameters(PairPotType1,PairPotType2,PotCoef)

implicit none
integer, parameter        :: PREC = 8
real,    parameter         :: ZERO = 0.0D0
character(1), intent(out) :: PairPotType1(10), PairPotType2(10)
real(PREC), intent(out)   :: PotCoef(16,10)
integer                   :: PPID, I, J
real(PREC)                :: R1, RCUT, R1SQ, POLY, SCL_R1, EXPTMP, DPOLY, DDPOLY
real(PREC)                :: DELTA, DELTA2, DELTA3, DELTA4

PotCoef = ZERO
PairPotType1(1) = 'H'
PairPotType2(1) = 'H'
PotCoef(1:10,1) = [38.512100000, 3.887860000, -37.769100000, 57.083500000, -34.512200000, &
  0.000000,0.000000,0.000000,0.900000,1.000000]

PairPotType1(2) = 'C'
PairPotType2(2) = 'C'
PotCoef(1:10,2) = [3.962931000,24.467772000,-51.156024000,39.031644000,-11.342979000, &
  0.000000,0.000000,0.000000,1.600000,1.700000] 

PairPotType1(3) = 'C'
PairPotType2(3) = 'H'
PotCoef(1:10,3) = [104.889589000,3.971095000,-23.823043000,26.408093000,-11.317522000, &
  0.000000,0.000000,0.000000,1.200000,1.300000]

PairPotType1(4) = 'N'
PairPotType2(4) = 'N'
PotCoef(1:10,4) = [43.228899000,15.004605000,-36.621777000,29.234888000,-8.912743000, &
  0.000000,0.000000,0.000000,1.600000,1.700000]

PairPotType1(5) = 'O'
PairPotType2(5) = 'O'
PotCoef(1:10,5) = [10.999870000,19.303033000,-45.747853000,37.946431000,-11.935755000, & 
  0.000000,0.000000,0.000000,1.500000,1.600000]

PairPotType1(6) = 'N'
PairPotType2(6) = 'C'
PotCoef(1:10,6) = [88.953762000,10.294988000,-27.706877000,22.101434000,-6.836438000, &
  0.000000,0.000000,0.000000,1.600000,1.700000]

PairPotType1(7) = 'N'
PairPotType2(7) = 'H'
PotCoef(1:10,7) = [0.625470000,28.081241000,-63.414297000,53.286361000,-17.352234000, &
  0.000000,0.000000,0.000000,1.300000,1.40000]

PairPotType1(8) = 'O'
PairPotType2(8) = 'N'
PotCoef(1:10,8) = [13.182426000,20.050322000,-46.806321000,38.206953000,-12.319656000, &
  0.000000,0.000000,0.000000,1.600000,1.700000]

PairPotType1(9) = 'O'
PairPotType2(9) = 'C'
PotCoef(1:10,9) = [0.944093000,30.116337000,-59.608215000,45.107654000,-13.178839000, &
  0.000000,0.000000,0.000000,1.600000,1.700000]

PairPotType1(10) = 'O'
PairPotType2(10) = 'H'
PotCoef(1:10,10) = [0.481176000,33.175383000,-81.158683000,74.935408000,-26.792315000, &
  0.000000,0.000000,0.000000,1.200000,1.300000]

!do I = 1,10
!do J = 1,10
!  PotCoef(I,J) = (floor(PotCoef(I,J)*1000.D0))/1000.D0
!enddo
!enddo



do PPID = 1,10
     R1 = PotCoef(9,PPID)
     RCUT = PotCoef(10,PPID)
     R1SQ = R1*R1

     POLY = R1*(PotCoef(2,PPID) + R1*(PotCoef(3,PPID) +  R1*(PotCoef(4,PPID) + R1*PotCoef(5,PPID))))
     SCL_R1 = PotCoef(1,PPID)*exp(POLY)
     EXPTMP = PotCoef(6,PPID)*exp(PotCoef(7,PPID)*(R1 - PotCoef(8,PPID)))
     PotCoef(11,PPID) = SCL_R1 + EXPTMP
     DPOLY = PotCoef(2,PPID) + 2.D0*PotCoef(3,PPID)*R1 + 3.D0*PotCoef(4,PPID)*R1SQ + 4.D0*PotCoef(5,PPID)*R1*R1SQ
     PotCoef(12,PPID) = DPOLY*SCL_R1 + PotCoef(7,PPID)*EXPTMP
     DDPOLY = 2.D0*PotCoef(3,PPID) + 6.D0*PotCoef(4,PPID)*R1 + 12.D0*PotCoef(5,PPID)*R1SQ
     PotCoef(13,PPID) = 0.5*((DPOLY*DPOLY + DDPOLY)*SCL_R1 + PotCoef(7,PPID)*PotCoef(7,PPID)*EXPTMP)

!    At the end of the join function:
     DELTA = RCUT - R1
     DELTA2 = DELTA*DELTA
     DELTA3 = DELTA2*DELTA
     DELTA4 = DELTA3*DELTA

     PotCoef(14,PPID) = (-1.D0/DELTA3)*(3.D0*PotCoef(13,PPID)*DELTA2 + 6.D0*PotCoef(12,PPID)*DELTA + 10.D0*PotCoef(11,PPID))
     PotCoef(15,PPID) = (1.D0/DELTA4)*(3.D0*PotCoef(13,PPID)*DELTA2 + 8.D0*PotCoef(12,PPID)*DELTA + 15.D0*PotCoef(11,PPID))

     PotCoef(16,PPID) = (-1.D0/(10.D0*DELTA3))*(6*PotCoef(15,PPID)*DELTA2 + 3.D0*PotCoef(14,PPID)*DELTA + PotCoef(13,PPID))
enddo

end subroutine LoadPairPotParameters
