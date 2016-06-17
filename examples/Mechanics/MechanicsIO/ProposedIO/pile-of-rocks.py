#!/usr/bin/env python

#
# Example of one object under gravity with one contactor and a ground
# using the Siconos proposed mechanics API
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics
import siconos.io.mechanics_io

import numpy
from scipy import constants

siconos.io.mechanics_io.use_proposed = True

options = siconos.mechanics.collision.bullet.BulletOptions()
options.worldScale = 1.0
options.breakingThreshold = 0.01

import pydoc
# Creation of the hdf5 file for input/output
with Hdf5() as io:

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.addNewtonImpactFrictionNSL('contact', mu=0.3)

    N = 1000
    sizes = numpy.linspace(0.7,3,N)
    w = int(numpy.sqrt(N)+1)
    s = 3
    x = numpy.linspace(-w*s,w*s,w)
    y = numpy.linspace(-w*s,w*s,w)
    z = numpy.linspace(0,w*s,w)
    numpy.random.shuffle(x)
    numpy.random.shuffle(y)
    numpy.random.shuffle(z)
    for i in range(N):
        pts = numpy.random.random((8,3))
        pts = (pts-pts.mean(0)) / (pts.max(0)-pts.min(0))
        io.addConvexShape('shape%d'%i, pts*sizes[i],
                          insideMargin=0.1, outsideMargin=0.0)
        x, y, z = numpy.random.random(3)*w*s
        io.addObject('obj%d'%i, [Contactor('shape%d'%i)],
                     translation=[x-w*s/2, y-w*s/2, z+2],
                     velocity=[0, 0, 0, 0, 0, 0],
                     mass=sizes[i]*0.1)

    # Definition of the ground shape
    io.addPrimitiveShape('Ground', 'Box', (1000, 1000, 1),
                         insideMargin=0.1, outsideMargin=0.0)

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.addObject('ground', [Contactor('Ground')],
                 translation=[0.0, 0.0, -0.5])

    # Definition of the shovel shape
    pts = numpy.array([(-10,10,0),
                       (-10,0,0),
                       (10,0,0),
                       (10,10,0),
                       (-10,0,1),
                       (10,0,1)])*10
    height = pts.mean(0)[2]
    io.addConvexShape('Shovel', pts-pts.mean(0),
                      insideMargin=0.01, outsideMargin=0.1)

    # # the shovel object made with the shovel shape
    io.addObject('shovel', [Contactor('Shovel')],
                 translation=[0.0, -150.0, height*1.0],
                 mass=100)

    io.addJoint('joint1', 'shovel',
                pivot_point=[0, 0, 0],
                axis=[0, 1, 0],
                joint_class='PrismaticJointR')

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with Hdf5(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.

    # print(pydoc.render_doc(io.run, "Help on %s"))

    from siconos.mechanics.collision import BodyDS, \
        BodyTimeStepping, SiconosContactor
    from siconos.mechanics.collision.bullet import BulletBroadphase

    def initialforces(body):
        weight = [0, 0, - body.scalarMass() * constants.g]
        body.setFExtPtr(weight)
        m = [0,0,0.0]
        body.setMExtPtr(m)

    from siconos.mechanics.collision.bullet import btQuaternion
    def myforces(name, body):
        if name=='shovel':
            weight = [0, 0, - body.scalarMass() * constants.g]
            q = body.q()
            #r = btQuaternion(q[4], q[5], q[6], q[3])
            #a = r.getAngle()
            #ax = r.getAxis()
            #print(r)
            #print(q[3:7])
            #print(a, ax.x(), ax.y(), ax.z())
            v = body.velocity()

            # PD controller
            weight[1] += (3.0-v[1])*1000
            #weight[2] += (40-q[2])*200 - v[2]*20
            body.setFExtPtr(weight)
            m = [0,0,0]
            body.setMExtPtr(m)

    io.run(with_timer=False,
           time_stepping=BodyTimeStepping,
           space_filter=lambda x: BulletBroadphase(x, options),
           body_class=BodyDS,
           shape_class=SiconosContactor,
           face_class=None,
           edge_class=None,
           t0=0,
           T=200,
           h=0.005,
           multipoints_iterations=True,
           theta=0.50001,
           Newton_max_iter=1,
           set_external_forces=initialforces,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=1000,
           tolerance=1e-4,
           numerics_verbose=False,
           body_callback=myforces,
           output_frequency=10)
