import sympy as sp
from sympy import Symbol
import numpy as np

q1 = np.array([Symbol('q10'), Symbol('q11'), Symbol('q12'), Symbol('q13')])
q2 = np.array([Symbol('q20'), Symbol('q21'), Symbol('q22'), Symbol('q23')])
cq2q10 = np.array([Symbol('_cq2q101'),Symbol('_cq2q102'),
                  Symbol('_cq2q103'),Symbol('_cq2q104')])
G1P0q1 = np.array([0, Symbol('_G1P0->getValue(0)'), Symbol('_G1P0->getValue(1)'),
                 Symbol('_G1P0->getValue(2)')])
G2P0q2 = np.array([0, Symbol('_G2P0->getValue(0)'), Symbol('_G2P0->getValue(1)'),
                 Symbol('_G2P0->getValue(2)')])

G1 = np.array([0, Symbol('X1'), Symbol('Y1'), Symbol('Z1')])
G2 = np.array([0, Symbol('X2'), Symbol('Y2'), Symbol('Z2')])

G1P1q1 = np.array([0, Symbol('_G1P1q1->getValue(0)'), Symbol('_G1P1q1->getValue(1)'),
                   Symbol('_G1P1q1->getValue(2)')])
G2P2q2 = np.array([0, Symbol('_G2P2q2->getValue(0)'), Symbol('_G2P2q2->getValue(1)'),
                   Symbol('_G2P2q2->getValue(2)')])

#TODO these might need to be rotated to q1 frame
P0P1q1 = np.array([0, Symbol('_axes[0]->getValue(0)'), Symbol('_axes[0]->getValue(1)'),
                   Symbol('_axes[0]->getValue(2)')])
P0P2q1 = np.array([0, Symbol('_axes[1]->getValue(0)'), Symbol('_axes[1]->getValue(1)'),
                   Symbol('_axes[1]->getValue(2)')])

# V00 = np.array([0, Symbol('_V[0][0]->getValue(0)'), Symbol('_V[0][0]->getValue(1)'),
#                Symbol('_V[0][0]->getValue(2)')])
# V01 = np.array([0, Symbol('_V[0][1]->getValue(0)'), Symbol('_V[0][1]->getValue(1)'),
#                Symbol('_V[0][1]->getValue(2)')])
# V10 = np.array([0, Symbol('_V[1][0]->getValue(0)'), Symbol('_V[1][0]->getValue(1)'),
#                Symbol('_V[1][0]->getValue(2)')])
# V11 = np.array([0, Symbol('_V[1][1]->getValue(0)'), Symbol('_V[1][1]->getValue(1)'),
#                Symbol('_V[1][1]->getValue(2)')])

qinv = lambda q: np.array([q[0],-q[1],-q[2],-q[3]])
qmul = lambda a,b: np.array([
         a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
         a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
         a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
         a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0]])

unrot = lambda V,q: qmul(qinv(q), qmul(V, q))
rot = lambda V,q: qmul(q, qmul(V, qinv(q)))

cross = lambda a,b: np.array([0]+list(np.cross(a[1:],b[1:])))

## Two-DS case

def getH(G1, q1, G2, q2):

    # World coordinates of all points
    P0q1 = G1 + rot(G1P0q1,q1)
    P0q2 = G2 + rot(G2P0q2,q2)
    # P1 = G1 + rot(G1P1q1,q1)
    # P2 = G2 + rot(G2P2q2,q2)

    # rot2to1 = qmul(q1,qinv(qmul(q2,cq2q10)))
    # P0q1 = G1P0q1
    # P0q2 = unrot((G2 + rot(G2P0q2,q2)) - G1, q1)

    # First part, position constraint (like knee joint, no free axis for now)
    HP = P0q1 - P0q2
    H1 = HP[1]
    H2 = HP[2]
    H3 = HP[3]

    # Second part, angle constraint, no rotation around the cross
    # product of P1P0 and P2P0', where P2P0' is P2P0 rotated around P1P0.
    rot2to1 = qmul(q1,qinv(qmul(q2,cq2q10)))

    # Measure rotation around first axis
    A = np.dot(P0P1q1, rot2to1) # = sin(angle)
    angle = 2*sp.atan2(rot2to1[0], A)
    return [H1,H2,H3,angle]

    # Use second axis rotated by that amount as basis for constrained
    # rotation vector: build a quaternion with P0P1 as axis, A as angle
    sA = sp.sin(angle/2)
    cA = sp.cos(angle/2)
    rA = np.hstack([[cA], sA*P0P1q1[1:]])

    V = rot(P0P2q1, rA)
    # H4 = np.dot(V, rot2to1)

    B = np.dot(V, rot2to1)
    angle1 = 2*sp.atan2(rot2to1[0], B)

    sB = sp.sin(angle/2)
    cB = sp.cos(angle/2)
    rB = np.hstack([[cB], sB*V[1:]])

    V = rot(rot(cross(P0P1q1,P0P2q1), rA),rB)
    H4 = np.dot(V, rot2to1)

    # V = unrot(cross(P0P1q1, P0P2q1), rot2to1)
    # H4 = np.dot(V, rot2to1)

    return [H1,H2,H3,H4]

H = getH(G1,q1,G2,q2)
dq = list(G1[1:])+list(q1)+list(G2[1:])+list(q2)
jachq = [[h.diff(d) for d in dq] for h in H]

with open('../pivot2_H.generated_c','w') as out:
    pre, exprs = sp.cse(H, optimizations = 'basic')
    for p in pre:
        print('  const double {} = {};'
              .format(p[0], str(sp.ccode(p[1]))),
              file=out)
    for n in range(len(H)):
        print('  y.setValue({}, {});'
              .format(n,str(sp.ccode(exprs[n]))),
              file=out)

with open('../pivot2_jachq_jd1d2.generated_c','w') as out:
    pre, exprs = sp.cse([x for y in jachq for x in y],
                        optimizations = 'basic')
    for p in pre:
        print('  const double {} = {};'
              .format(p[0], str(sp.ccode(p[1]))),
              file=out)
    n = 0
    for i in range(len(H)):
        for j in range(len(dq)):
            print('  _jachq->setValue({}, {}, {});'.format(
                i, j, str(sp.ccode(exprs[n]))), file=out)
            n += 1

## One-DS case

H = getH(G1,q1,np.array([0,0,0,0]),np.array([1,0,0,0]))
dq = list(G1[1:])+list(q1)
jachq = [[h.diff(d) for d in dq] for h in H]

with open('../pivot2_jachq_jd1.generated_c','w') as out:
    pre, exprs = sp.cse([x for y in jachq for x in y],
                        optimizations = 'basic')
    for p in pre:
        print('  const double {} = {};'
              .format(p[0], str(sp.ccode(p[1]))),
              file=out)
    n = 0
    for i in range(len(H)):
        for j in range(len(dq)):
            print('  _jachq->setValue({}, {}, {});'.format(
                i, j, str(sp.ccode(exprs[n]))), file=out)
            n += 1

import os
os.system('touch ../Pivot2JointR.cpp')
