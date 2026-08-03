"""Fisher-discriminant empirical mold: the optimal projective base mold.
w = Sigma^-1 (mu1 - mu0), means from label-free GMM clustering on PCA10,
Sigma = regularized noise covariance. Score = (patch - mid) . w."""
import csv, os
import numpy as np
from sklearn.mixture import GaussianMixture

SIDE = 17
COORDS = np.arange(-0.32, 0.321, 0.04)
TT, UU = np.meshgrid(COORDS, COORDS, indexing='ij')
DISK = (TT**2 + UU**2) <= (0.32**2)
GRID = np.arange(SIDE * SIDE).reshape(SIDE, SIDE)
DISK_IDX = GRID[DISK]

def load_patches(path, prefix):
    out = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            key = (os.path.basename(r["file"]), int(r["lobe"]))
            vals, ok = [], True
            for i in range(1, 290):
                v = r.get(f"{prefix}_p{i:03d}", "NA")
                try: vals.append(float(v))
                except ValueError: ok = False; break
            out[key] = np.array(vals).reshape(SIDE, SIDE) if ok else None
    return out

def flip_u_disk(flat):
    g = np.zeros(SIDE * SIDE)
    g[DISK_IDX] = flat
    return g.reshape(SIDE, SIDE)[:, ::-1].ravel()[DISK_IDX]

def fit_fisher(tr_X, tr_amps):
    """GMM on PCA10 -> means -> Fisher vector in patch space."""
    Xc = tr_X - tr_X.mean(0)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    V10 = Vt[:10]
    Z = Xc @ V10.T
    gm = GaussianMixture(n_components=2, covariance_type="full", random_state=0,
                         reg_covar=1e-3, max_iter=300).fit(Z)
    a = gm.predict(Z)
    mu_p = np.array([Xc[a == c].mean(0) if (a == c).any() else Xc.mean(0) for c in range(2)])
    ma = [tr_amps[a == c].mean() for c in range(2)]
    order = np.argsort(ma)  # low amp = GlcN
    g0, g1 = mu_p[order[0]], mu_p[order[1]]
    # Fisher vector in the 10-dim latent space
    cov = np.cov(Z.T) + 1e-2 * np.eye(10)
    z0 = (g0 - Xc.mean(0)) @ V10.T
    z1 = (g1 - Xc.mean(0)) @ V10.T
    w_z = np.linalg.solve(cov, z1 - z0)
    w_p = V10.T @ w_z  # back to patch space (centered)
    mid = Xc.mean(0)
    return w_p, mid, (g0, g1)

def score(x, w_p, mid):
    return float((x - mid) @ w_p)

def main():
    import sys
    if len(sys.argv) < 4:
        print("Usage: empirical_fisher_mold.py PATCHES_TSV PREFIX OUT_PREFIX")
        print("Writes OUT_PREFIX.tsv (full-data Fisher) and OUT_PREFIX_cv.tsv")
        print("(half-split cross-validated Fisher margins). Label-free: patch")
        print("clustering + amplitude mapping only.")
        sys.exit(1)
    path, prefix, out_prefix = sys.argv[1], sys.argv[2], sys.argv[3]
    patches = load_patches(path, prefix)
    keys = sorted(k for k in patches if patches[k] is not None)
    X = np.array([patches[k][DISK] for k in keys])
    amps = np.array([patches[k][8, 8] for k in keys])

    def write_scores(tag, w_p, mid, out):
        with open(out, "w") as g:
            g.write("file\tlobe\tpredicted\tconfidence\tscore\n")
            for key, x in zip(keys, X):
                s = max(score(x, w_p, mid), score(flip_u_disk(x), w_p, mid))
                g.write(f"{key[0]}\t{key[1]}\t{1 if s < 0 else 0}\t{abs(s):.8f}\t{s:.6f}\n")
        print(f"{tag}: {out}")

    w_p, mid, _ = fit_fisher(X, amps)
    write_scores("full", w_p, mid, out_prefix + ".tsv")

    even = [i for i in range(len(keys)) if keys[i][1] % 2 == 0]
    odd = [i for i in range(len(keys)) if keys[i][1] % 2 == 1]
    w_e, mid_e, _ = fit_fisher(X[even], amps[even])
    w_o, mid_o, _ = fit_fisher(X[odd], amps[odd])
    with open(out_prefix + "_cv.tsv", "w") as g:
        g.write("file\tlobe\tpredicted\tconfidence\tscore\n")
        for i, key in enumerate(keys):
            w, m = (w_o, mid_o) if i in even else (w_e, mid_e)
            s = max(score(X[i], w, m), score(flip_u_disk(X[i]), w, m))
            g.write(f"{key[0]}\t{key[1]}\t{1 if s < 0 else 0}\t{abs(s):.8f}\t{s:.6f}\n")
    print(f"cv: {out_prefix}_cv.tsv")


if __name__ == "__main__":
    main()
