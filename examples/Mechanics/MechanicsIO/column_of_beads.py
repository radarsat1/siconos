#!/usr/bin/env python

#
# Example of two cubes, one with a convex shape, one with a primitive
# shape.
#

from siconos.mechanics.contact_detection.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics

n_beads = 1000

object_type = 'ConvexHull' # 'Sphere', 'Box'

# Creation of the hdf5 file for input/output
with Hdf5() as io:
    for n in range(n_beads):
        # Definition of a bead as a sphere (diameter 1.0)
        io.addPrimitiveShape('Sphere%d'%n, 'Sphere', [0.5])

        # Definition of a "bead" as a box (dimensions 1.0)
        io.addPrimitiveShape('Box%d'%n, 'Box', [1.0]*3)

        # Definition of a cube as a convex shape (dimensions 1.0)
        io.addConvexShape('Hull%d'%n,
                          [(-0.5, 0.5, -0.5), (-0.5, -0.5, -0.5),
                           (-0.5, -0.5, 0.5), (-0.5, 0.5, 0.5),
                           (0.5, 0.5, 0.5), (0.5, 0.5, -0.5),
                           (0.5, -0.5, -0.5), (0.5, -0.5, 0.5)])

    # The bead object made with an unique Contactor : the sphere shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0
    for n in range(n_beads):
        io.addObject('sphere%d'%n, [Contactor('Sphere%d'%n)],
                     translation=[0, 0, 1.5*(n+1)],
                     mass=1)
        io.addObject('box%d'%n, [Contactor('Box%d'%n)],
                     translation=[-2, 0, 1.5*(n+1)],
                     mass=1)
        io.addObject('hull%d'%n, [Contactor('Hull%d'%n)],
                     translation=[2, 0, 1.5*(n+1)],
                     mass=1)

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.addPrimitiveShape('Ground', 'Box', (12, 4, 1))
    io.addObject('ground', [Contactor('Ground')],
                 translation=[0, 0, 0])

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.addNewtonImpactFrictionNSL('contact', mu = 1.0, e = 0.0)

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.

nstep=10000
step=0.001
with Hdf5(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(with_timer=False,
           time_stepping=None,
           space_filter=None,
           body_class=None,
           shape_class=None,
           face_class=None,
           edge_class=None,
           gravity_scale=1,
           t0=0,
           T=nstep*step,
           h=step,
           multipoints_iterations=True,
           theta=0.50001,
           Newton_max_iter=1,
           set_external_forces=None,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=10000,
           tolerance=1e-4,
           numerics_verbose=False,
           output_frequency=1)
