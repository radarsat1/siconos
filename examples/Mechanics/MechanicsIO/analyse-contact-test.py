#!/usr/bin/env python3

import h5py, json, pandas
from collections import defaultdict
from matplotlib import pyplot as plt
import numpy as np

h = h5py.File('contact-test-data-run2.hdf5', mode='r')

df_z_min = pandas.DataFrame()
df_z_std = pandas.DataFrame()
df_contacts_mean = pandas.DataFrame()
df_contacts_std = pandas.DataFrame()
df_time = pandas.DataFrame()

stats = defaultdict(Stats)
for g in h:
    config   = json.loads(h[g].attrs['config'])
    contacts = h[g]['contacts']
    simtime  = h[g]['simtime']
    walltime = h[g]['walltime']
    z        = h[g]['z']
    rng      = np.argmax(np.array(simtime)>3)

    t = config['threshold']
    df_z_min = df_z_min.append( {t: np.min(z)}, ignore_index=True )
    df_z_std = df_z_std.append( {t: np.std(z)}, ignore_index=True )
    df_contacts_mean = df_contacts_mean.append({t: np.mean(contacts)}, ignore_index=True )
    df_contacts_std = df_contacts_std.append({t: np.std(contacts)}, ignore_index=True )
    df_time = df_time.append( {t: np.mean(walltime[-1] - walltime[rng])}, ignore_index=True )

#     plt.plot(simtime[rng:],np.array(contacts[rng:]) - config['size']/2)
# plt.show()

df_z_min.sort_index(axis='columns', inplace=True)
df_z_std.sort_index(axis='columns', inplace=True)
df_contacts_mean.sort_index(axis='columns', inplace=True)
df_contacts_std.sort_index(axis='columns', inplace=True)
df_time.sort_index(axis='columns', inplace=True)

plt.subplot(3,2,1)
plt.semilogx(df_contacts_mean.columns, df_contacts_mean.mean())
plt.title('contacts (mean)')
plt.subplot(3,2,2)
plt.semilogx(df_contacts_std.columns, df_contacts_std.mean())
plt.title('contacts (std)')
plt.subplot(3,2,3)
plt.semilogx(df_z_min.columns, df_z_min.mean())
plt.title('z (min)')
plt.subplot(3,2,4)
plt.semilogx(df_z_std.columns, df_z_std.mean())
plt.title('z (std)')
plt.subplot(3,2,5)
plt.semilogx(df_time.columns, df_time.mean())
plt.title('time')
plt.subplot(3,2,6)

def normed(df):
    return (df.mean() - df.mean().mean()) / df.mean().std()

plt.semilogx(df_contacts_mean.columns,
             - (normed(df_contacts_mean))
             + (normed(df_contacts_std))
             - (normed(df_z_min))
             + (normed(df_z_std))
             + (normed(df_time)))
plt.title('combined')
plt.show()
