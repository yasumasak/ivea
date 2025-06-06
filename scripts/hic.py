import numpy as np
import scipy.sparse as ssp
import pandas as pd
import time, os

def get_hic_file(chromosome, hic_dir, allow_vc=True, hic_type="juicebox"):
    if hic_type == "juicebox":
        hic_file = os.path.join(hic_dir, chromosome, chromosome + ".KRobserved.gz")
        hic_norm = os.path.join(hic_dir, chromosome, chromosome + ".KRnorm.gz")

        is_vc = False
        if allow_vc and not hic_exists(hic_file):
            hic_file = os.path.join(hic_dir, chromosome, chromosome + ".VCobserved.gz")
            hic_norm = os.path.join(hic_dir, chromosome, chromosome + ".VCnorm.gz")

            if not hic_exists(hic_file):
                RuntimeError("Could not find KR or VC normalized hic files")
            else:
                print("  Could not find KR normalized hic file. Using VC normalized hic file")
                is_vc = True

        print("  Using: " + hic_file)
        return hic_file, hic_norm, is_vc
    elif hic_type == "bedpe":
        hic_file = os.path.join(hic_dir, chromosome, chromosome + ".bedpe.gz")

        return hic_file, None, None


def hic_exists(file):
    if not os.path.exists(file):
        return False
    elif file.endswith('gz'):
        # gzip file still have some size. This is a hack
        return (os.path.getsize(file) > 100)
    else:
        return (os.path.getsize(file) > 0)


def load_hic(hic_file, hic_norm_file, hic_is_vc, hic_type, hic_resolution, tss_hic_contribution, window, min_window, gamma, interpolate_nan=True, apply_diagonal_bin_correction=True):
    print("  Loading HiC data")

    if hic_type == 'juicebox':
        HiC_sparse_mat = hic_to_sparse(hic_file, hic_norm_file, hic_resolution)
        HiC = process_hic(hic_mat = HiC_sparse_mat,
                            hic_norm_file = hic_norm_file,
                            hic_is_vc = hic_is_vc,
                            resolution = hic_resolution,
                            tss_hic_contribution = tss_hic_contribution,
                            window = window,
                            min_window = min_window,
                            gamma = gamma,
                            interpolate_nan = interpolate_nan,
                            apply_diagonal_bin_correction = apply_diagonal_bin_correction)

    elif hic_type == 'bedpe':
        HiC = pd.read_csv(hic_file, sep="\t", names = ['chr1','x1','x2','chr2','y1','y2','name','hic_contact'])

    return HiC


def hic_to_sparse(filename, norm_file, resolution, hic_is_doubly_stochastic=False):

    hic = pd.read_table(filename, names=["bin1", "bin2", "hic_contact"],
                        header=None, engine='c', memory_map=True)

    # Verify our assumptions
    assert np.all(hic.bin1 <= hic.bin2)

    # Need load norms here to know the dimensions of the hic matrix
    norms = pd.read_csv(norm_file, header=None)
    hic_size = norms.shape[0]

    # Convert to sparse matrix in CSR (compressed sparse row) format, chopping
    # down to HiC bin size. note that conversion to scipy sparse matrices
    # accumulates repeated indices, so this will do the right thing.
    row = np.floor(hic.bin1.values / resolution).astype(int)
    col = np.floor(hic.bin2.values / resolution).astype(int)
    dat = hic.hic_contact.values

    #JN: Need both triangles in order to compute row/column sums to make double stochastic.
    #If juicebox is upgraded to return DS matrices, then can remove one triangle
    #TO DO: Remove one triangle when juicebox is updated.
    # we want a symmetric matrix.  Easiest to do that during creation, but have to be careful of diagonal
    if not hic_is_doubly_stochastic:
        mask = (row != col)  # off-diagonal
        row2 = col[mask]  # note the row/col swap
        col2 = row[mask]
        dat2 = dat[mask]

        # concat and create
        row = np.hstack((row, row2))
        col = np.hstack((col, col2))
        dat = np.hstack((dat, dat2))

    return ssp.csr_matrix((dat, (row, col)), (hic_size, hic_size))


def process_hic(hic_mat, hic_norm_file, hic_is_vc, resolution, tss_hic_contribution, window, min_window=0, hic_is_doubly_stochastic=False, apply_diagonal_bin_correction=True, interpolate_nan=True, gamma=None, kr_cutoff = .25):
    # Make doubly stochastic.
    # Juicer produces a matrix with constant row/column sums. But sum is not 1 and is variable across chromosomes

    if not hic_is_doubly_stochastic and not hic_is_vc:
        #Any row with Nan in it will sum to nan
        #So need to calculate sum excluding nan
        temp = hic_mat
        temp.data = np.nan_to_num(temp.data, copy=False)
        sums = temp.sum(axis = 0)
        sums = sums[~np.isnan(sums)]

        # For KR-normalized matrix
        assert(np.max(sums[sums > 0])/np.min(sums[sums > 0]) < 1.001)
        mean_sum = np.mean(sums[sums > 0])
        if abs(mean_sum - 1) < .001:
            print('  HiC matrix has row sums of {}, continuing without making doubly stochastic'.format(mean_sum))
        else:
            print('  HiC matrix has row sums of {}, making doubly stochastic...'.format(mean_sum))
            hic_mat = hic_mat.multiply(1/mean_sum)

    # Adjust diagonal of matrix based on neighboring bins
    # First and last rows need to be treated differently
    if apply_diagonal_bin_correction:
        last_idx = hic_mat.shape[0] - 1
        nonzero_diag = hic_mat.nonzero()[0][hic_mat.nonzero()[0] == hic_mat.nonzero()[1]]
        nonzero_diag = list(set(nonzero_diag) - set(np.array([last_idx])) - set(np.array([0])))

        for ii in nonzero_diag:
            hic_mat[ii,ii] = max(hic_mat[ii,ii-1], hic_mat[ii,ii+1]) * tss_hic_contribution / 100

        if hic_mat[0,0] != 0:
            hic_mat[0, 0] = hic_mat[0,1] * tss_hic_contribution / 100

        if hic_mat[last_idx, last_idx] != 0:
            hic_mat[last_idx, last_idx] = hic_mat[last_idx, last_idx - 1] * tss_hic_contribution / 100

    # Any entries with low KR norm entries get set to NaN. These will be interpolated below
    hic_mat = apply_kr_threshold(hic_mat, hic_norm_file, kr_cutoff)

    # Remove lower triangle
    if not hic_is_vc:
        hic_mat = ssp.triu(hic_mat)
    else:
        hic_mat = process_vc(hic_mat)

    # Turn into dataframe
    hic_mat = hic_mat.tocoo(copy=False)
    hic_df = pd.DataFrame({'bin1': hic_mat.row, 'bin2': hic_mat.col, 'hic_contact': hic_mat.data})

    # Prune to window
    hic_df = hic_df.loc[np.logical_and(abs(hic_df['bin1'] - hic_df['bin2']) <= window/resolution, abs(hic_df['bin1'] - hic_df['bin2']) >= min_window/resolution)]
    print("  HiC has {} rows after windowing between {} and {}".format(hic_df.shape[0], min_window, window))

    # Fill NaN
    # NaN in the KR normalized matrix are not zeros. They are entries where the KR algorithm did not converge (or low KR norm)
    # So need to fill these. Use powerlaw.
    # Not ideal obviously but the scipy interpolation algos are either very slow or don't work since the nan structure implies that not all nans are interpolated
    if interpolate_nan:
        nan_loc = np.isnan(hic_df['hic_contact'])
        hic_df.loc[nan_loc,'hic_contact'] = get_powerlaw_at_distance(abs(hic_df.loc[nan_loc,'bin1'] - hic_df.loc[nan_loc,'bin2']) * resolution, gamma)

    return(hic_df)


def apply_kr_threshold(hic_mat, hic_norm_file, kr_cutoff):
    # Convert all entries in the HiC matrix corresponding to low kr norm entries to NaN
    # Note that in scipy sparse matrix multiplication 0*nan = 0
    # So this doesn't convert 0's to nan only nonzero to nan

    norms = np.loadtxt(hic_norm_file)
    norms[norms < kr_cutoff] = np.nan
    norms[norms >= kr_cutoff] = 1
    norm_mat = ssp.dia_matrix(( 1.0/norms, [0]), (len(norms), len(norms)))

    return norm_mat * hic_mat * norm_mat


def get_powerlaw_at_distance(distances, gamma, min_distance=5000, scale=None):
    assert(gamma > 0)

    #The powerlaw is computed for distances > 5kb. We don't know what the contact freq looks like at < 5kb.
    #So just assume that everything at < 5kb is equal to 5kb.
    #TO DO: get more accurate powerlaw at < 5kb
    distances = np.clip(distances, min_distance, np.Inf)
    log_dists = np.log(distances + 1)

    #Determine scale parameter
    #A powerlaw distribution has two parameters: the exponent and the minimum domain value 
    #In our case, the minimum domain value is always constant (equal to 1 HiC bin) so there should only be 1 parameter
    #The current fitting approach does a linear regression in log-log space which produces both a slope (gamma) and a intercept (scale)
    #Empirically there is a linear relationship between these parameters (which makes sense since we expect only a single parameter distribution)
    #It should be possible to analytically solve for scale using gamma. But this doesn't quite work since the hic data does not actually follow a power-law
    #So could pass in the scale parameter explicity here. Or just kludge it as I'm doing now
    #TO DO: Eventually the pseudocount should be replaced with a more appropriate smoothing procedure.

    #4.80 and 11.63 come from a linear regression of scale on gamma across 20 hic cell types at 5kb resolution. Do the params change across resolutions?
    if scale is None:
        scale = -4.80 + 11.63 * gamma

    powerlaw_contact = np.exp(scale + -1*gamma * log_dists)

    return(powerlaw_contact)


def process_vc(hic):
    # For a vc normalized matrix, need to make rows sum to 1.
    # Assume rows correspond to genes and cols to enhancers

    row_sums = hic.sum(axis = 0)
    row_sums[row_sums == 0] = 1
    norm_mat = ssp.dia_matrix((1.0 / row_sums, [0]), (row_sums.shape[1], row_sums.shape[1]))

    # Left multiply to operate on rows
    hic = norm_mat * hic

    return(hic)

