
from __future__ import print_function
import os,sys
import numpy as np
import math

np.set_printoptions(precision=3)

from siconos.mechanics.collision.tools import Contactor
from siconos.mechanics.joints import cast_PivotJointR
from siconos.io.mechanics_io import Hdf5
from siconos.kernel import SiconosVector, BlockVector, changeFrameAbsToBody
import siconos.kernel as Kernel

# Use of Pivot2 with and without first axis free (suspension)

# Creation of the hdf5 file for input/output
with Hdf5() as io:

    ## Definition of a bar and a wheel
    io.addPrimitiveShape('Point', 'Sphere', (0.01,))
    point = [Contactor('Point')]
    io.addPrimitiveShape('Bar', 'Box', (0.5, 0.1, 0.1))
    io.addPrimitiveShape('VBar', 'Box', (0.1, 0.1, 0.5))
    io.addPrimitiveShape('Coupler', 'Box', (0.01, 0.03, 0.01))
    io.addPrimitiveShape('Handle', 'Box', (0.05, 0.05, 0.05))
    io.addPrimitiveShape('Wheel', 'Cylinder', (0.3, 0.1))
    bar = [Contactor('Bar')]
    vbar = [Contactor('VBar')]
    wheel = [#Contactor('Wheel'),
             #Contactor('Handle', relative_translation=[0,0.075,0.275]),
             #Contactor('Handle', relative_translation=[0,0.075,-0.275]),
             Contactor('Coupler', relative_translation=[0.01,-0.5,0],
                       relative_orientation=[(0,0,1),np.pi/2]),
             Contactor('Bar', relative_translation=[0,-0.2,0],
                       relative_orientation=[(0,0,1),np.pi/2])]
    coupler = [Contactor('Coupler', relative_translation=[0,0.01,0]),
               Contactor('Coupler', relative_translation=[0,0,0.01],
                         relative_orientation=[(1,0,0),np.pi/2])]

    #vel = [1,0,0,1,3,0]
    vel = [1,0,0,0,3,0]

    # io.addObject('point1', point, [-0.5,0,0.5])
    io.addObject('point2', point, [-0.5,-0.5,0.5])

    ## Definition of the ground
    io.addPrimitiveShape('Ground', 'Box', (5, 7, 0.1))
    io.addObject('ground', [Contactor('Ground')], [0,0,-0.05])
    # io.addNewtonImpactFrictionNSL('contact', mu=0.3, e=0.0)

    ## Fixed object (two ds), one pivot2
    io.addObject('bar1', vbar, [-0.5,-0.5,0.8], mass=1)#orientation=[(0,1,0),np.pi/2], mass=1)
    io.addObject('wheel1', wheel, [-0.5,0,0.5], mass=1, velocity=vel)
    io.addJoint('joint1', 'bar1', None, None, None, 'FixedJointR')
    io.addJoint('joint2','bar1','wheel1',[-0.5,-0.5,0.5],[[0,0,1],[0,1,0]],'Pivot2JointR')

    ## Fixed object (two ds), two pivots with coupler
    io.addObject('bar2', vbar, [0.5,-0.5,0.8], mass=1)#orientation=[(0,1,0),np.pi/2], mass=1)
    io.addObject('wheel2', wheel, [0.5,0,0.5], mass=1, velocity=vel)
    io.addObject('coupler1', coupler, [0.5,-0.5,0.5], mass=1)
    io.addJoint('joint3', 'bar2', None, None, None, 'FixedJointR')
    io.addJoint('joint4','bar2','coupler1',[0.5,-0.5,0.5],[0,0,1],'PivotJointR')
    io.addJoint('joint5','coupler1','wheel2',[0.5,-0.5,0.5],[0,1,0],'PivotJointR')

# Load and run the simulation
with Hdf5(mode='r+') as io:
    io.run(t0=0,
           T=5,
           h=0.01,
           theta=0.5,
           Newton_max_iter=1,
           itermax=1000,
           tolerance=1e-12,
           projection_itermax=3,
           projection_tolerance=1e-5,
           projection_tolerance_unilateral=1e-5,
           time_stepping=Kernel.TimeSteppingDirectProjection,
           osi=Kernel.MoreauJeanDirectProjectionOSI,
           set_external_forces = lambda x: None,
           )
