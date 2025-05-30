# Alon Diament, Tuller Lab, June 2022.

import numpy as np
from scipy.io import savemat, loadmat


def save_SA(path, SA):
    SA.pop('mask', None)
    SA.pop('homologs', None)
    np.savez_compressed(path, **SA)


def load_SA(path):
    with np.load(path) as npz:
        SA = {k: v for k, v in npz.items()}
    return SA


def save_matlab_SA(path, SA):
    """ generating MATLAB variables: SA (n x 6 matrix), ref (string cell).
    """
    savemat(path,
            {'SA': np.vstack([SA['pos'] + 1, SA['ind'] + 1,
                              np.ones(SA['pos'].size),
                              SA['pos_from_stop']]).T,
             'ref': np.array([[np.array(r)] for r in SA['ref']], dtype=object)})


def load_matlab_SA(path):
    """ expected MATLAB variables: SA (n x 6 matrix), ref (string cell).
    """
    SA_mat = loadmat(path)
    return {'ref': [r[0] for r in np.ravel(SA_mat['ref'])],
            'ind': SA_mat['SA'][:, 1] - 1,
            'pos': SA_mat['SA'][:, 0] - 1,
            'pos_from_stop': SA_mat['SA'][:, 3]}
