#!/usr/bin/env python3

# Run with:
#   blender -P blender_siconos.py

# Tested on Blender 2.79b

import bpy
import bmesh
import math
import h5py
import numpy as np

data = h5py.File('cube.hdf5', mode='r')['data']

instances = {}

# Clear default cube
default_cube = bpy.data.objects.get("Cube")
if default_cube is not None:
    bpy.data.objects.remove(default_cube)

# Load scene
def load_scene():
    for oname in data['input']:
        for instname in data['input'][oname]:
            o = data['input'][oname]
            inst = o[instname]
            shapename = inst.attrs['shape_name']
            shape = data['ref'][shapename]
            mesh = None
            if shape.attrs['type'] == 'convex':
                verts = shape[:]
                mesh = bpy.data.meshes.new(instname+'_'+shapename)
                mesh.from_pydata(verts,[],[])
                bm = bmesh.new()
                bm.from_mesh(mesh)
                bmesh.ops.convex_hull(bm, input=bm.verts)
                bm.to_mesh(mesh)
            elif shape.attrs['type'] == 'primitive':
                if shape.attrs['primitive'] == 'Box':
                    w,h,d = np.transpose(shape)
                    verts = [(-w/2,-h/2,-d/2),(-w/2,h/2,-d/2),(w/2,h/2,-d/2),(w/2,-h/2,-d/2),
                             (-w/2,-h/2,d/2),(-w/2,h/2,d/2),(w/2,h/2,d/2),(w/2,-h/2,d/2)]
                    faces = [(0,1,2,3), (4,5,6,7), (0,4,5,1), (1,5,6,2), (2,6,7,3), (3,7,4,0)]
                    mesh = bpy.data.meshes.new(instname+':'+shapename)
                    mesh.from_pydata(verts,[],faces)
            if mesh is not None:
                obj = bpy.data.objects.new(instname, mesh)
                obj.location = o.attrs['translation'][:]
                obj.rotation_mode = 'QUATERNION'
                mesh.update(calc_edges=True)
                bpy.context.scene.objects.link(obj)
                if o.attrs['id'] not in instances:
                    instances[o.attrs['id']] = []
                instances[o.attrs['id']].append(obj)

def position_static_objects():
    for p in data['static']:
        for k in range(p.shape[0]):
            i = int(p[1])
            for obj in instances[i]:
                # TODO shape offsets
                obj.location = tuple(p[2:5])
                obj.rotation_quaternion = tuple(p[5:9])

def position_dynamic_objects(time):
    pos = data['dynamic']
    k = np.searchsorted(pos[:,0], time)
    while pos[k,0] <= time:
        i = int(pos[k,1])
        for obj in instances[i]:
            # TODO shape offsets
            obj.location = tuple(pos[k,2:5])
            obj.rotation_quaternion = tuple(pos[k,5:9])
        # TODO invisible objects that were not touched
        k += 1

# every frame change, this function is called.
def frame_handler(scene):
    time = scene.frame_current / scene.render.fps
    position_dynamic_objects(time)

def setup():
    load_scene()
    position_static_objects()
    position_dynamic_objects(0.0)
    bpy.app.handlers.frame_change_pre.append(frame_handler)

setup()
