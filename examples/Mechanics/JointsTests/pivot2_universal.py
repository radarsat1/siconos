
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

# A demonstration that Pivot2 with orthogonal axes, Universal, and two
# sequential Pivot joints all behave the same.

# Creation of the hdf5 file for input/output
with Hdf5() as io:

    # Definition of a bar and a coupler
    io.addPrimitiveShape('Bar', 'Box', (1, 0.1, 0.1))
    io.addPrimitiveShape('Coupler', 'Box', (0.1, 0.01, 0.01))
    coupler = [Contactor('Coupler',relative_orientation=[(0,0,1),np.pi/2]),
               Contactor('Coupler',relative_orientation=[(0,1,0),np.pi/2])]

    # Definition of the ground
    io.addPrimitiveShape('Ground', 'Box', (5, 7, 0.1))
    io.addObject('ground', [Contactor('Ground')], [0,0,-0.05])

    # Function to express some repeated cases
    vel = [0,0,0,3,3,10]
    obj, jnt = 1, 1
    def twojoints(y, args):
        global obj, jnt

        # Two test cases:
        # background (one ds)
        # io.addObject('bar%d'%obj, [Contactor('Bar')], [1.2,y,2], mass=1, velocity=vel)
        # io.addJoint('joint%d'%jnt, 'bar%d'%obj, None, *args(1.8,y))
        # obj += 1
        # jnt += 1

        # Fixed object (two ds)
        io.addObject('bar%d'%obj, [Contactor('Bar')], [0,y,2], mass=1)
        io.addObject('bar%d'%(obj+1), [Contactor('Bar')], [-1.2,y,2], mass=1, velocity=vel)
        io.addJoint('joint%d'%jnt, 'bar%d'%obj, None, None, None, 'FixedJointR')
        io.addJoint('joint%d'%(jnt+1), 'bar%d'%obj, 'bar%d'%(obj+1), *args(-0.6,y))
        obj += 2
        jnt += 2

    ## Regular pivot (hinge) joint.
    # twojoints(2, lambda x,y: ([x,y,2], [0,1,0], 'PivotJointR'))

    ## The knee (ball) joint rotates in one more direction compared to
    ## the pivot-based joints.
    # twojoints(1, lambda x,y: ([x,y,2], [], 'KneeJointR'))

    ## The pivot2 joint: like two pivot joints but no need for an extra
    ## simulation body in-between.  We make the two axes orthogonal to
    ## simulate the universal joint.
    twojoints(0, lambda x,y: ([x,y,2], [[0,1,0], [0,0,1]], 'Pivot2JointR'))

    ## The universal joint: like two pivot joints with orthogonal axes,
    ## but no need for an extra simulation body in-between.
    # twojoints(-1, lambda x,y: ([[x,y,2]], [[0,1,0], [0,0,1]], 'UniversalJointR'))

    ## Two sequential pivot joints. we make the two axes orthogonal to
    ## simulate the universal joint.
    # background (one ds)
    jt = 'PivotJointR'
    y = -2
    # io.addObject('bar%d'%obj, [Contactor('Bar')], [1.2,y,2], mass=1, velocity=vel)
    # io.addObject('coupler%d'%obj, coupler, [1.8,y,2], mass=1, velocity=vel)
    # io.addJoint('joint%d'%jnt, 'coupler%d'%obj, None, [1.8,y,2], [0,1,0], jt)
    # io.addJoint('joint%d'%(jnt+1), 'bar%d'%obj, 'coupler%d'%obj, [1.8,y,2], [0,0,1], jt)
    # obj += 1
    # jnt += 2

    # Fixed object (two ds)
    io.addObject('bar%d'%obj, [Contactor('Bar')], [0,y,2], mass=1)
    io.addObject('coupler%d'%obj, coupler, [-0.6,y,2], mass=1)
    io.addObject('bar%d'%(obj+1), [Contactor('Bar')], [-1.2,y,2], mass=1, velocity=vel)
    io.addJoint('joint%d'%jnt, 'bar%d'%obj, None, None, None, 'FixedJointR')
    io.addJoint('joint%d'%(jnt+1), 'bar%d'%obj, 'coupler%d'%obj, [-0.6,y,2], [0,1,0], jt)
    io.addJoint('joint%d'%(jnt+2), 'coupler%d'%obj, 'bar%d'%(obj+1), [-0.6,y,2], [0,0,1], jt)
    obj += 2
    jnt += 3

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
