#!/usr/bin/env python

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics
import siconos.io.mechanics_io
from siconos.mechanics.collision.convexhull import ConvexHull
import numpy as np
import time, json, h5py, multiprocessing

parameters = {'size': [10, 1.0, 0.1, 0.01],
              'margin': [0.01, 0.1, 0.4],
              'groundmargin': [0.01, 0.1, 1.0],
              'height': [1, 2, 3, 20],
              'angle': [0, 0.01, 0.1],
              'shape': ['tetra', 'cube', 'box'],
              'timestep': [0.01, 0.001],
              'perturbation': [3,4,5],
              'newton': [1, 2, 5],
              'moving': [True,False],
              'boththresholds': [True,False],
              'threshold': [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3,
                            0.01, 0.02, 0.03, 0.04, 0.1, 0.2, 0.3]}

# Note: Why does this configuration crash with normq > 0 assertion?
# (NewtonEulerDS.cpp:335)

# config = {'angle': 0, 'newton': 5, 'threshold': 1e-10,
# 'perturbation': 4, 'height': 3, 'groundmargin': 1.0,
# 'boththresholds': False, 'shape': 'tetra', 'size': 0.1, 'timestep':
# 0.01, 'margin': 0.01, 'moving': True}

def generate_configuration():
    valid = False
    while not valid:
        config = {}
        for k in parameters.keys():
            config[k] = parameters[k][np.random.randint(len(parameters[k]))]
        valid = (config['margin'] < (config['size']/2)
                 and config['height'] > (config['size']/2))
    return config

def design_world(config):
    with Hdf5() as io:
        velocity = [0,0,0,0,0,0]
        if config['moving']:
            velocity=[-4.121433, -0.220321, 1.259349,
                      -0.124028, 0.231024, 0.710144]

        a = np.random.uniform(-np.pi, np.pi)
        orientation = [(np.cos(a), np.sin(a), 0), float(config['angle'])]
        print('=== orientation',orientation)

        if config['shape'] is 'tetra':
            pts = np.array([(-1.0, 1.0, 1.0),
                            (1.0, -1.0, 1.0),
                            (-1.0, -1.0, 1.0),
                            (0.0, 0.0, -1.0)])
            pts *= config['size']
            io.addConvexShape('Tetra', pts - pts.mean(0),
                              insideMargin=config['margin'], outsideMargin=0.0)

            ch = ConvexHull(pts)
            inertia,volume=ch.inertia(ch.centroid())

            io.addObject('obj', [Contactor('Tetra')],
                         translation=[0, 0, config['height']],
                         orientation=orientation,
                         velocity=velocity,
                         mass=1, inertia=inertia)

        elif config['shape'] is 'cube':
            pts = np.array([
                (-1.0, 1.0, -1.0),
                (-1.0, -1.0, -1.0),
                (-1.0, -1.0, 1.0),
                (-1.0, 1.0, 1.0),
                (1.0, 1.0, 1.0),
                (1.0, 1.0, -1.0),
                (1.0, -1.0, -1.0),
                (1.0, -1.0, 1.0)])
            pts *= config['size']/2
            io.addConvexShape('Cube', pts - pts.mean(0),
                              insideMargin=config['margin'], outsideMargin=0.0)

            ch = ConvexHull(pts)
            inertia,volume=ch.inertia(ch.centroid())

            io.addObject('obj', [Contactor('Cube')],
                         translation=[0, 0, config['height']],
                         orientation=orientation,
                         velocity=velocity,
                         mass=1, inertia=inertia)

        elif config['shape'] is 'box':
            io.addPrimitiveShape('Cube', 'Box', np.array([1,1,1])*config['size'],
                                 insideMargin=config['margin'], outsideMargin=0.0)

            io.addObject('obj', [Contactor('Cube')],
                         translation=[0, 0, config['height']],
                         orientation=orientation,
                         velocity=velocity,
                         mass=1)

        io.addPrimitiveShape('Ground', 'Box', (100, 100, 100),
                             insideMargin=config['groundmargin'], outsideMargin=0.0)

        io.addNewtonImpactFrictionNSL('contact', mu=0.01, e=0.01)

        io.addObject('ground', [Contactor('Ground')],
                     translation=[0, 0, -50-config['size']/2])


def run_world(config):
    class controller():
        def initialize(self,io):
            self.io = io
            topo = io._model.nonSmoothDynamicalSystem().topology()
            self.obj = topo.getDynamicalSystem('obj')
            self.trace = []
            self.t0 = time.time()
        def step(self):
            tm = self.io.currentTime()
            contacts = (self.io._model.simulation()
                        .oneStepNSProblem(0).getSizeOutput()//3)
            self.trace.append([time.time()-self.t0, tm, contacts, self.obj.q()[2]])

    ctrl = controller()

    options = siconos.mechanics.collision.bullet.SiconosBulletOptions()
    options.worldScale = 1.0
    options.perturbationIterations = config['perturbation']
    options.minimumPointsPerturbationThreshold = config['perturbation']
    options.manifoldContactBreakingThreshold = config['threshold']
    if config['boththresholds']:
        options.contactBreakingThreshold = config['threshold']
    with Hdf5(mode='r+') as io:
        io.run(options=options,
               controller=ctrl,
               t0=0,
               T=5,
               h=config['timestep'],
               theta=0.50001,
               Newton_max_iter=config['newton'],
               solver=Numerics.SICONOS_FRICTION_3D_NSGS,
               itermax=100000,
               tolerance=1e-8,
               output_frequency=None)
    return ctrl.trace

def save_config_trace(h, counter, config, trace):
    g = h.create_group('%04d'%counter)
    g.attrs['config'] = json.dumps(config)
    g.create_dataset('walltime',  data=[t[0] for t in trace])
    g.create_dataset('simtime',   data=[t[1] for t in trace])
    g.create_dataset('contacts',  data=[t[2] for t in trace])
    g.create_dataset('z',         data=[t[3] for t in trace])

def one_run(counter, config):
    def task(counter, config):
        print('===',config)
        design_world(config)
        trace = run_world(config)
        h = h5py.File('contact-test-data.hdf5', mode='a')
        save_config_trace(h, counter, config, trace)
        print('task done')
    p = multiprocessing.Process(target=task, args=(counter,config))
    p.start()
    return p

if __name__=="__main__":
    for counter in range(10000):
        print()
        print('=== {} ==='.format(counter))
        config = generate_configuration()
        p = one_run(counter, config)
        p.join()
        print('process joined')
