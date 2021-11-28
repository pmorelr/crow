from cylinderModule import *

dt = 0.01
tend = 10.

cylinder = LIAproblem2D(100,"./Results/Re100/")
cylinder.nonlinearTimestepper_Euler(dt, tend)
