"""Compare octree-based and gold-standard eigenvector reconstruction.

Port of ``mainEigReconstructMethods.m``.

Usage::

    python main_eig_reconstruct_methods.py                       # default file list
    python main_eig_reconstruct_methods.py a.nrrd b.nrrd c.nrrd  # custom volumes
    python main_eig_reconstruct_methods.py --plot --seed 0 ...

Output NRRD files are written to the current directory (or ``--out``).
"""
import argparse

from gold_eig_reconstruct import test_gold_eig_reconstruct
from octree_eig_reconstruct import test_octree_eig_reconst_interp

# toy examples
# FILE_LIST = ['./Toothuc0.20.nrrd', './Toothucderiv10.20.nrrd', './Toothucderiv20.20.nrrd']
FILE_LIST = ['./Pf25cpUs0.20.nrrd', './TCf25cpUs0.20.nrrd',
             './QVAPORf25cpUs0.20.nrrd', './PRECIPf25cpUs0.20.nrrd']

# full-sized
# FILE_LIST = ['../data/Toothuc.nrrd', '../data/Toothucderiv1.nrrd', '../data/Toothucderiv2.nrrd']
# FILE_LIST = ['../data/Pf25cpUs.nhdr', '../data/TCf25cpUs.nhdr',
#              '../data/QVAPORf25cpUs.nhdr', '../data/PRECIPf25cpUs.nhdr']


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('files', nargs='*', default=FILE_LIST, help='input .nrrd/.nhdr volumes')
    parser.add_argument('--val-diff', type=float, default=0.02)
    parser.add_argument('--out', default='.', help='output directory')
    parser.add_argument('--seed', type=int, default=None, help='random seed for Monte-Carlo PCA')
    parser.add_argument('--plot', action='store_true', help='show intermediate figures')
    parser.add_argument('--skip-gold', action='store_true', help='skip the gold standard')
    args = parser.parse_args()

    # octree eigenvectors computation
    pca_oct_vols = test_octree_eig_reconst_interp(args.files, args.val_diff,
                                                  show_plots=args.plot, output_dir=args.out,
                                                  seed=args.seed)
    # gold-standard
    if not args.skip_gold:
        pca_gold_vols = test_gold_eig_reconstruct(args.files, output_dir=args.out,
                                                  seed=args.seed)

    if args.plot:
        import matplotlib.pyplot as plt
        plt.show()


if __name__ == '__main__':
    main()
