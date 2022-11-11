!pip3 install uproot awkward
import uproot, scipy, math, os, random, time, collections
import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns
import awkward as ak
from tqdm import tqdm
import pandas as pd

f_fluk= uproot.open('/content/drive/MyDrive/data/alignment/fluka_ke.root')
f_coll= uproot.open('/content/drive/MyDrive/data/alignment/kfalignment_data_iter0.root')

collision = f_coll['trackParam'].arrays(['fitParam_x',
                                         'fitParam_y',
                                         'fitParam_z',
                                         'fitParam_px',
                                         'fitParam_py',
                                         'fitParam_pz',
                                         'fitParam_chi2',
                                         'fitParam_charge',
                                         'fitParam_align_local_residual_x_sp',
                                         'fitParam_align_global_residual_y_sp'], library='ak')

fluka = f_fluk['trackParam'].arrays(['fitParam_x',
                                     'fitParam_y',
                                     'fitParam_z',
                                     'fitParam_px',
                                     'fitParam_py',
                                     'fitParam_pz',
                                     'fitParam_chi2',
                                     'fitParam_charge',
                                     'fitParam_align_local_residual_x_sp',
                                     'fitParam_align_global_residual_y_sp'], library='ak')

plt.figure(figsize=(14, 14))
plt.scatter(ak.flatten(collision['fitParam_x']), ak.flatten(collision['fitParam_y']), marker='.', s=.05, label='Collision')
plt.scatter(ak.flatten(fluka['fitParam_x']), ak.flatten(fluka['fitParam_y']), marker='.', s=.2, label='Fluka')
plt.xlim(-150, 150)
plt.ylim(-150, 150)
plt.xlabel('x (mm)', fontsize=16)
plt.ylabel('y (mm)', fontsize=16)
plt.legend(fontsize=14)
plt.title('x-y scatter plot', fontsize=20)
plt.savefig('xy.png')

plt.figure(figsize=(13, 9))
bin = np.lin
plt.hist(ak.flatten(collision['fitParam_z']), bins=101, alpha=0.5, density=True, label='Collision')
plt.hist(ak.flatten(fluka['fitParam_z']), bins=101, alpha=0.5, density=True, label='Fluka')
plt.xlabel('z (mm)', fontsize=18)
plt.ylabel('Number of Events (normalized)', fontsize=18)
plt.legend(fontsize=14)
plt.title('Z Position', fontsize=22)
plt.savefig('z.png')

# ty
plt.figure(figsize=(25, 10))
plt.subplot(1, 2, 1)
bin = np.linspace(-0.1, 0.1, 101)
fluka_tx, coll_tx = ak.flatten(fluka['fitParam_px'])/ak.flatten(fluka['fitParam_pz']), ak.flatten(collision['fitParam_px'])/ak.flatten(collision['fitParam_pz'])
plt.hist(fluka_tx, bins=bin, alpha=0.5, density=True, label=f'Fluka (MC) tx µ: {round(np.mean(fluka_tx), 6)}, σ: {round(np.std(fluka_tx), 6)}')
plt.hist(coll_tx, bins=bin, alpha=0.5, density=True, label=f'Collision tx µ: {round(np.mean(coll_tx), 6)}, σ: {round(np.std(coll_tx), 6)}')

plt.xlabel('tx', fontsize=18)
plt.ylabel('Number of Events (normalized)', fontsize=18)
plt.title('Collision tx', fontsize=22)
plt.xticks(fontsize=15)
plt.legend(fontsize=18)

# ty
plt.subplot(1, 2, 2)
fluka_ty, coll_ty = ak.flatten(fluka['fitParam_py'])/ak.flatten(fluka['fitParam_pz']), ak.flatten(collision['fitParam_py'])/ak.flatten(collision['fitParam_pz'])
plt.hist(fluka_ty, bins=bin, alpha=0.5, density=True, label=f'Fluka (MC) ty µ: {round(np.mean(fluka_ty), 6)}, σ: {round(np.std(fluka_ty), 6)}')
plt.hist(coll_ty, bins=bin, alpha=0.5, density=True, label=f'Collision ty µ: {round(np.mean(coll_ty), 6)}, σ: {round(np.std(coll_tx), 6)}')

plt.xlabel('ty', fontsize=18)
plt.ylabel('Number of Events (normalized)', fontsize=18)
plt.title('Collsion ty', fontsize=22)
plt.xticks(fontsize=15)

plt.legend(fontsize=18)
plt.savefig('txty.png', dpi=fig.dpi)

fluka_p = np.sqrt(ak.flatten(fluka['fitParam_px'])**2 + ak.flatten(fluka['fitParam_py'])**2 + ak.flatten(fluka['fitParam_pz'])**2)
coll_p = np.sqrt(ak.flatten(collision['fitParam_px'])**2 + ak.flatten(collision['fitParam_py'])**2 + ak.flatten(collision['fitParam_pz'])**2)
plt.figure(figsize=(13, 9))
plt.hist(fluka_p, bins=np.linspace(0, 2e3, 201), alpha=0.5, density=True, label='Fluka')
plt.hist(coll_p, bins=np.linspace(0, 2e3, 201), alpha=0.5, density=True, label='Collision')
plt.xlabel('p (GeV)', fontsize=18)
plt.ylabel('Number of Events (normalized)', fontsize=18)
plt.legend(fontsize=14)
plt.title('Momentum', fontsize=22)
plt.savefig('momentum.png', dpi=fig.dpi)

plt.figure(figsize=(13, 9))
plt.hist(ak.flatten(fluka['fitParam_chi2']), bins=np.linspace(0, 150, 151), alpha=0.5, density=True, label='Fluka')
plt.hist(ak.flatten(collision['fitParam_chi2']), bins=np.linspace(0, 150, 151), alpha=0.5, density=True, label='Collision')
plt.xlabel('Chi2', fontsize=18)
plt.ylabel('Number of Events (normalized)', fontsize=18)
plt.legend(fontsize=14)
plt.title('Chi2', fontsize=22)
plt.savefig('chi2.png', dpi=fig.dpi)
