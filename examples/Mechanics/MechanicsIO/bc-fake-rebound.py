#!/usr/bin/env python

#
# Example of a sphere bouncing but with damped velocity generated
# "artificially" using boundary conditions.
#

from siconos.mechanics.collision.tools import Contactor
import siconos.numerics as Numerics
import siconos.kernel as Kernel
from siconos.io import mechanics_io
import numpy as np

# Creation of the hdf5 file for input/output
with mechanics_io.Hdf5() as io:

    # Definition of a sphere
    io.addPrimitiveShape('Sphere', 'Sphere', (2,),
                         insideMargin=0.2, outsideMargin=0.3)

    # Definition of the ground shape
    io.addPrimitiveShape('Ground', 'Box', (10, 10, 0.1),
                         insideMargin=0.05, outsideMargin=0.1)

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.addNewtonImpactFrictionNSL('contact', mu=0.1, e=0.9)

    # The sphere object made with an unique Contactor : the sphere shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0
    io.addObject('sphere', [Contactor('Sphere')],
                 translation=[0, 0, 4],
                 velocity=[0, 0, 0, 0, 0, 0],
                 mass=1)

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.addObject('ground', [Contactor('Ground')],
                 translation=[0, 0, -0.1])

    io.addBoundaryCondition('myBC', 'sphere', indices=[],
                            bc_class='DampedBounceBC')

# Example: "Fake" rebound calculations in the boundary condition
@mechanics_io.namespace
class DampedBounceBC(Kernel.BoundaryCondition):
    def __init__(self, params=None):
        # Initialize the class with a vector of the maximum length we need.
        # There are 3 indices (0,1,2) for linear velocity.
        super(self.__class__, self).__init__([0,1,2,3,4,5],[0,0,0,0,0,0])

        # State
        self.collision = False
        self.velocity = 10.0

    def updateBoundaryConditionsForDS(self, ds, vfree):
        # We are in collision if norm(p(1)) is non-zero.
        print('p',ds.p(1))        # Sum of reaction forces computed by solver
        print('v',ds.velocity())  # Incoming velocity 6-vector
        print('q',ds.q())         # Incoming position 7-vector
        if np.linalg.norm(ds.p(1)) > 0.001:
            self._velocityIndices = [0,1,2,3,4,5]
            self._prescribedVelocity.zero()
            self._prescribedVelocity.setValue(0, 0.4)
            if not self.collision:
                # Use max(v,0) to ensure we don't "pull down" on the ball
                self.velocity = np.max([-ds.velocity()[2] * 0.95, 0])
            self._prescribedVelocity.setValue(2, self.velocity)
            # We use a flag to keep track of whether it is the first
            # timestep of the bounce.
            self.collision = True
        else:
            # Don't affect any velocity parts by default.
            if self.collision:
                self.velocity *= 0.9
            self.collision = False
            self._velocityIndices = []

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with mechanics_io.Hdf5(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    io.run(t0=0,
           T=20,
           h=0.001,
           theta=0.50001,
           Newton_max_iter=1,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=1000,
           tolerance=1e-5)
