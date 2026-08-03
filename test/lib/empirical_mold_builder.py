"""Empirical molds: cluster the experimental patches into 2 shape centroids
(label-free), map by the physical amplitude convention, score by NCC."""
import csv, os, collections
import numpy as np

SIDE = 17
COORDS = np.arange(-0.32, 0.321, 0.04)
TT, UU = np.meshgrid(COORDS, COORDS, indexing='ij')
DISK = (TT**2 + UU**2) <= (0.32**2)

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

def kmeans2(X, seeds=(0, 1, 2, 3, 4), iters=100):
    best = None
    for seed in seeds:
        rng = np.random.default_rng(seed)
        idx = rng.choice(len(X), 2, replace=False)
        mu = X[idx].copy()
        for _ in range(iters):
            d = ((X[:, None, :] - mu[None, :, :]) ** 2).sum(-1)
            a = d.argmin(1)
            new = np.array([X[a == k].mean(0) if (a == k).any() else mu[k] for k in range(2)])
            if np.abs(new - mu).max() < 1e-9:
                mu = new
                break
            mu = new
        d = ((X[:, None, :] - mu[None, :, :]) ** 2).sum(-1)
        inertia = d.min(1).sum()
        if best is None or inertia < best[0]:
            best = (inertia, mu, d.argmin(1))
    return best[1], best[2]

def main():
    import sys
    if len(sys.argv) < 4:
        print("Usage: empirical_mold_builder.py PATCHES_TSV PREFIX OUT_TSV")
        sys.exit(1)
    path, prefix, out = sys.argv[1], sys.argv[2], sys.argv[3]
    patches = load_patches(path, prefix)
    keys = [k for k in patches if patches[k] is not None]
    X = np.array([patches[k][DISK] for k in keys])
    print(f"patches: {len(keys)}, dims: {X.shape[1]}")

    mu, assign = kmeans2(X)
    amps = np.array([patches[k][8, 8] for k in keys])
    mean_amp = [amps[assign == c].mean() for c in range(2)]
    print(f"amplitude moyenne par cluster: {mean_amp}")
    glcn = mu[int(np.argmin(mean_amp))]
    glcnac = mu[int(np.argmax(mean_amp))]

    def ncc(a, b):
        av = a - a.mean(); bv = b - b.mean()
        d = np.sqrt((av * av).sum() * (bv * bv).sum())
        return float((av * bv).sum() / d) if d > 0 else 0.0

    with open(out, "w") as g:
        g.write("file\tlobe\tpredicted\tconfidence\tncc0\tncc1\n")
        for key in sorted(keys):
            p = patches[key][DISK]
            n0 = ncc(p, glcn)
            n1 = ncc(p, glcnac)
            margin = n1 - n0
            g.write(f"{key[0]}\t{key[1]}\t{1 if margin < 0 else 0}\t{-margin:.8f}\t{n0:.6f}\t{n1:.6f}\n")
    print(f"empirical molds scored -> {out}")
    return (glcn, glcnac)


if __name__ == "__main__":
    main()
