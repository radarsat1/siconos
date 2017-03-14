#!/usr/bin/env python

#
# Example of one object under gravity with one contactor and a ground
# using the Siconos proposed mechanics API
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io import mechanics_io
import siconos.numerics as Numerics
import siconos.kernel as Kernel
import numpy as np
from siconos.mechanics.collision.convexhull import ConvexHull

# Creation of the hdf5 file for input/output
with mechanics_io.Hdf5() as io:

    # Definition of a tetrahedron as a convex shape
    import numpy
    pts = numpy.array([(-1.0, 1.0, 1.0),
                       (1.0, -1.0, 1.0),
                       (-1.0, -1.0, 1.0),
                       (0.0, 0.0, -1.0)])
    io.addConvexShape('Tetra', pts - pts.mean(0),
                      insideMargin=0.01, outsideMargin=0.0)

    # Definition of the ground shape
    io.addPrimitiveShape('Ground', 'Box', (10, 10, 1),
                         insideMargin=0.01, outsideMargin=0.0)

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.addNewtonImpactFrictionNSL('contact', mu=0.01, e=0.7,
                                  collision_group1=1,
                                  collision_group2=2)

    # computation of inertia and volume
    ch = ConvexHull(pts)
    inertia,volume=ch.inertia(ch.centroid())

    # The tetra object made with an unique Contactor : the tetrahedron
    # shape.  As a mass is given, it is a dynamic system involved in
    # contact detection and in the simulation.  With no group id
    # specified the Contactor belongs to group 0
    io.addObject('tetra', [Contactor('Tetra', collision_group=1)],
                 translation=[0, 0, 4],
                 velocity=[0, 0, 0, 0, 0, 0],
                 mass=1, inertia=inertia)

    io.addBoundaryCondition('myBC', 'tetra', indices=[],
                            bc_class='FirstCollisionNoRotateBC')

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.addObject('ground', [Contactor('Ground', collision_group=2)],
                 translation=[0, 0, -0.1])

# Example: Suppress the angular velocity on first bounce
@mechanics_io.namespace
class FirstCollisionNoRotateBC(Kernel.BoundaryCondition):
    def __init__(self, params=None):
        # Initialize the class with a vector of the maximum length we need.
        # There are 3 indices (3,4,5) for angular velocity.
        super(self.__class__, self).__init__([3,4,5],[0,0,0])

        # State
        self.collision = False
        self.count = 0

    def updateBoundaryConditionsForDS(self, ds, vfree):
        # We are in collision if norm(p(1)) is non-zero.
        # Increment count on first timestep of a collision.
        c = self.count
        if np.linalg.norm(ds.p(1)) > 0.001:
            if not self.collision:
                self.count += 1
            self.collision = True
        else:
            self.collision = False

        # If this is the first bounce, specify the angular velocity
        # indices to set corresponding parts of the velocity to zero.
        if self.collision and self.count == 1:
            self._velocityIndices = [3,4,5]
        else:
            # Otherwise, don't replace any velocity parts.
            self._velocityIndices = []

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with mechanics_io.Hdf5(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.

    io.run(with_timer=False,
           t0=0,
           T=20,
           h=0.005,
           theta=0.50001,
           Newton_max_iter=1,
           set_external_forces=None,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=100000,
           tolerance=1e-8,
           numerics_verbose=False,
           output_frequency=None)
