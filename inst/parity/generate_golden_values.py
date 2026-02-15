#!/usr/bin/env python3
"""Generate golden reference values from PyGSP for rgsp parity testing.

Run:  cd inst/parity && python generate_golden_values.py

Outputs JSON files to ../extdata/parity/
"""

import json
import os
import sys

import numpy as np
from scipy import sparse

# Ensure the local pygsp is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'pygsp'))

import pygsp
from pygsp import graphs, filters, learning, reduction

OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'extdata', 'parity')
os.makedirs(OUT_DIR, exist_ok=True)

np.random.seed(42)


def to_list(x):
    """Convert numpy arrays to plain Python lists for JSON serialization."""
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, (np.floating, np.integer)):
        return float(x)
    return x


def save_golden(filename, data):
    """Save golden values as JSON."""
    data['_meta'] = {
        'pygsp_version': pygsp.__version__,
        'numpy_version': np.__version__,
        'schema_version': 1,
        'generator': 'generate_golden_values.py',
    }
    path = os.path.join(OUT_DIR, filename)
    with open(path, 'w') as f:
        json.dump(data, f, indent=2)
    print(f"  Wrote {path}")


# =============================================================================
# A. Graphs: eigenvalues and lmax for deterministic graphs
# =============================================================================
print("A. Generating golden_graphs.json ...")

graph_specs = {
    'Ring_10':         lambda: graphs.Ring(10),
    'Path_8':          lambda: graphs.Path(8),
    'Grid2d_3_4':      lambda: graphs.Grid2d(3, 4),
    'FullConnected_6': lambda: graphs.FullConnected(6),
    'Torus_4_5':       lambda: graphs.Torus(4, 5),
}

graph_data = {}
for name, constructor in graph_specs.items():
    g = constructor()
    g.compute_fourier_basis()
    g.estimate_lmax()
    graph_data[name] = {
        'n': int(g.N),
        'eigenvalues': to_list(np.sort(g.e)),
        'lmax': float(g.lmax),
        'lmax_exact': float(np.max(g.e)),
    }

save_golden('golden_graphs.json', graph_data)


# =============================================================================
# B. GFT: Ring(8) -- gft coefficients (abs sorted for sign ambiguity)
# =============================================================================
print("B. Generating golden_gft.json ...")

g = graphs.Ring(8)
g.compute_fourier_basis()

np.random.seed(42)
x_gft = np.random.randn(8)
gft_coeffs = g.gft(x_gft)
# Round-trip check
x_rec = g.igft(gft_coeffs)

gft_data = {
    'graph': 'Ring_8',
    'signal': to_list(x_gft),
    'gft_abs_sorted': to_list(np.sort(np.abs(gft_coeffs))),
    'round_trip_error': float(np.max(np.abs(x_rec - x_gft))),
}

save_golden('golden_gft.json', gft_data)


# =============================================================================
# C. Heat filter: kernel at eigenvalues, filtered output, lmax
# =============================================================================
print("C. Generating golden_filters.json ...")

filter_data = {}

# C1: Ring(12), scale=0.5
g = graphs.Ring(12)
g.compute_fourier_basis()
g.estimate_lmax()
scale = 0.5
filt = filters.Heat(g, scale=scale)

np.random.seed(42)
x_heat = np.random.randn(12)
y_heat = filt.filter(x_heat)

# Kernel evaluated at eigenvalues: exp(-scale * e / lmax)
kernel_at_eig = np.exp(-scale * g.e / g.lmax)

filter_data['heat_ring12'] = {
    'graph': 'Ring_12',
    'scale': scale,
    'lmax': float(g.lmax),
    'eigenvalues': to_list(np.sort(g.e)),
    'kernel_at_eigenvalues': to_list(kernel_at_eig),
    'signal': to_list(x_heat),
    'filtered_output': to_list(y_heat.flatten()),
    'note': 'PyGSP kernel: exp(-scale * lambda / lmax). rgsp uses exp(-t * lambda), so t = scale / lmax.',
}

# C2: Grid2d(4,4), scale=1.0
g = graphs.Grid2d(4, 4)
g.compute_fourier_basis()
g.estimate_lmax()
scale = 1.0
filt = filters.Heat(g, scale=scale)

np.random.seed(42)
x_heat2 = np.random.randn(16)
y_heat2 = filt.filter(x_heat2)

kernel_at_eig2 = np.exp(-scale * g.e / g.lmax)

filter_data['heat_grid2d_4_4'] = {
    'graph': 'Grid2d_4_4',
    'scale': scale,
    'lmax': float(g.lmax),
    'eigenvalues': to_list(np.sort(g.e)),
    'kernel_at_eigenvalues': to_list(kernel_at_eig2),
    'signal': to_list(x_heat2),
    'filtered_output': to_list(y_heat2.flatten()),
}

save_golden('golden_filters.json', filter_data)


# =============================================================================
# D. MexicanHat bank: Ring(20), Nf=4
# =============================================================================
print("D. Generating golden_filterbanks.json ...")

filterbank_data = {}

g = graphs.Ring(20)
g.compute_fourier_basis()
g.estimate_lmax()
Nf = 4
filt = filters.MexicanHat(g, Nf=Nf)

# Kernel evaluations at eigenvalues (per band)
kernel_evals = filt.evaluate(g.e)

filterbank_data['mexicanhat_ring20'] = {
    'graph': 'Ring_20',
    'Nf': Nf,
    'lmax': float(g.lmax),
    'eigenvalues': to_list(np.sort(g.e)),
    'n_bands': len(kernel_evals),
    'scales': to_list(np.array(filt.scales)),
    'kernel_evals': [to_list(k) for k in kernel_evals],
    'note': 'PyGSP MexicanHat bandpass: scale*x*exp(-scale*x). Lowpass: 1.2*exp(-1)*exp(-(x/0.4/lmin)^4).',
}


# =============================================================================
# E. Meyer bank: Path(16), Nf=3
# =============================================================================
g = graphs.Path(16)
g.compute_fourier_basis()
g.estimate_lmax()
Nf = 3
filt = filters.Meyer(g, Nf=Nf)

kernel_evals_meyer = filt.evaluate(g.e)

filterbank_data['meyer_path16'] = {
    'graph': 'Path_16',
    'Nf': Nf,
    'lmax': float(g.lmax),
    'eigenvalues': to_list(np.sort(g.e)),
    'n_bands': len(kernel_evals_meyer),
    'scales': to_list(np.array(filt.scales)),
    'kernel_evals': [to_list(k) for k in kernel_evals_meyer],
}

save_golden('golden_filterbanks.json', filterbank_data)


# =============================================================================
# F. Wave kernel: Ring(12) -- for documentation/property tests
# =============================================================================
print("F. Generating golden_wave.json ...")

g = graphs.Ring(12)
g.compute_fourier_basis()
g.estimate_lmax()
filt_wave = filters.Wave(g, time=5, speed=1.0)

wave_kernel_evals = filt_wave.evaluate(g.e)

wave_data = {
    'graph': 'Ring_12',
    'time': 5,
    'speed': 1.0,
    'lmax': float(g.lmax),
    'eigenvalues': to_list(np.sort(g.e)),
    'kernel_at_eigenvalues': [to_list(k) for k in wave_kernel_evals],
    'note': 'PyGSP wave: cos(t*arccos(1-speed^2*x/lmax/2)). rgsp: cos(t*sqrt(x)). Different formulas.',
}

save_golden('golden_wave.json', wave_data)


# =============================================================================
# G. Chebyshev: Ring(12) Heat(0.5) m=30
# =============================================================================
print("G. Generating golden_chebyshev.json ...")

g = graphs.Ring(12)
g.compute_fourier_basis()
g.estimate_lmax()
scale = 0.5
filt = filters.Heat(g, scale=scale)

np.random.seed(42)
x_cheby = np.random.randn(12)

# Filter with Chebyshev approximation (method='chebyshev', order=30)
y_cheby = filt.filter(x_cheby, method='chebyshev', order=30)
# Filter with exact method for reference
y_exact = filt.filter(x_cheby, method='exact')

cheby_data = {
    'graph': 'Ring_12',
    'scale': scale,
    'lmax': float(g.lmax),
    'signal': to_list(x_cheby),
    'filtered_chebyshev_m30': to_list(y_cheby.flatten()),
    'filtered_exact': to_list(y_exact.flatten()),
    'chebyshev_order': 30,
}

save_golden('golden_chebyshev.json', cheby_data)


# =============================================================================
# H. Differential operators: Ring(10), Grid2d(3,3)
# =============================================================================
print("H. Generating golden_differential.json ...")

diff_data = {}

# H1: Ring(10)
g = graphs.Ring(10)
g.compute_differential_operator()
D = g.D

np.random.seed(42)
x_diff = np.random.randn(10)
grad_out = g.grad(x_diff)

diff_data['ring10'] = {
    'graph': 'Ring_10',
    'n_nodes': int(g.N),
    'D_shape_pygsp': list(D.shape),  # PyGSP: N x E convention
    'n_edges': int(max(D.shape)),    # edges = larger dimension
    'signal': to_list(x_diff),
    'grad_output': to_list(grad_out.flatten()),
    'note_D': 'PyGSP D is N x E. rgsp D is E x N (transposed convention).',
}

# H2: Grid2d(3,3)
g = graphs.Grid2d(3, 3)
g.compute_differential_operator()

np.random.seed(42)
x_diff2 = np.random.randn(9)
grad_out2 = g.grad(x_diff2)
div_grad_out = g.div(grad_out2)
L = g.L
Lx = L.dot(x_diff2)

diff_data['grid2d_3_3'] = {
    'graph': 'Grid2d_3_3',
    'n_nodes': int(g.N),
    'D_shape_pygsp': list(g.D.shape),  # PyGSP: N x E convention
    'n_edges': int(max(g.D.shape)),    # edges = larger dimension
    'signal': to_list(x_diff2),
    'grad_output': to_list(grad_out2.flatten()),
    'div_grad_output': to_list(div_grad_out.flatten()),
    'Lx': to_list(Lx),
    'note': 'In PyGSP: div(grad(x)) = L*x. In rgsp: div(grad(x)) = -L*x (sign convention).',
}

save_golden('golden_differential.json', diff_data)


# =============================================================================
# I. Tikhonov: Ring(20), tau=5.0
# =============================================================================
print("I. Generating golden_learning.json ...")

g = graphs.Ring(20)
g.compute_fourier_basis()
g.compute_laplacian()

np.random.seed(42)
x_tik = np.random.randn(20)
tau = 5.0

# PyGSP regression_tikhonov with mask=all-ones
mask = np.ones(g.N, dtype=bool)
y_tik = learning.regression_tikhonov(g, x_tik, mask, tau)

# Laplacian quadratic form
Lq_input = float(x_tik @ g.L.toarray() @ x_tik)
Lq_output = float(y_tik @ g.L.toarray() @ y_tik)

learning_data = {
    'tikhonov_ring20': {
        'graph': 'Ring_20',
        'tau': tau,
        'signal': to_list(x_tik),
        'output': to_list(y_tik),
        'laplacian_qf_input': Lq_input,
        'laplacian_qf_output': Lq_output,
        'note': 'PyGSP: (M + tau*L)x = M*y with M=I. rgsp: (I + tau*L)x = y. Equivalent when mask=all-ones.',
    }
}

save_golden('golden_learning.json', learning_data)


# =============================================================================
# J. Kron reduction: Ring(12), keep=[0,2,4,6,8,10]
# =============================================================================
print("J. Generating golden_reduction.json ...")

g = graphs.Ring(12)
g.compute_fourier_basis()
g.compute_laplacian()

keep_idx = [0, 2, 4, 6, 8, 10]  # 0-indexed for Python
g_red = reduction.kron_reduction(g, keep_idx)

# Compute eigenvalues of reduced graph
if hasattr(g_red, 'compute_fourier_basis'):
    g_red.compute_fourier_basis()
    red_eigs = np.sort(g_red.e)
else:
    L_red = g_red if sparse.issparse(g_red) else g_red.L
    red_eigs = np.sort(np.linalg.eigvalsh(L_red.toarray() if sparse.issparse(L_red) else L_red))

reduction_data = {
    'kron_ring12': {
        'graph': 'Ring_12',
        'keep_0indexed': keep_idx,
        'keep_1indexed': [k + 1 for k in keep_idx],
        'reduced_n': len(keep_idx),
        'reduced_eigenvalues': to_list(red_eigs),
        'note': 'PyGSP uses 0-indexed keep. rgsp uses 1-indexed keep.',
    }
}

save_golden('golden_reduction.json', reduction_data)


# =============================================================================
# K. Frame bounds: Ring(10) Heat(1.0), Ring(16) Regular
# =============================================================================
print("K. Generating golden_frames.json ...")

frames_data = {}

# K1: Heat frame bounds
g = graphs.Ring(10)
g.compute_fourier_basis()
g.estimate_lmax()

filt = filters.Heat(g, scale=1.0)
bounds = filt.estimate_frame_bounds()

frames_data['heat_ring10'] = {
    'graph': 'Ring_10',
    'filter': 'Heat(scale=1.0)',
    'A': float(bounds[0]),
    'B': float(bounds[1]),
    'lmax': float(g.lmax),
}

# K2: Regular frame bounds
g = graphs.Ring(16)
g.compute_fourier_basis()
g.estimate_lmax()

filt = filters.Regular(g)
bounds = filt.estimate_frame_bounds()

# Also get kernel evaluations
kernel_evals_reg = filt.evaluate(g.e)

frames_data['regular_ring16'] = {
    'graph': 'Ring_16',
    'filter': 'Regular',
    'A': float(bounds[0]),
    'B': float(bounds[1]),
    'lmax': float(g.lmax),
    'n_filters': len(kernel_evals_reg),
    'kernel_evals': [to_list(k) for k in kernel_evals_reg],
}

save_golden('golden_frames.json', frames_data)


# =============================================================================
# L. Tight frames: Ring(16) -- Regular, Held, Simoncelli, Papadakis
# =============================================================================
print("L. Generating tight frame kernel evals ...")

g = graphs.Ring(16)
g.compute_fourier_basis()
g.estimate_lmax()

tight_data = {}
for name, FiltClass in [('Regular', filters.Regular),
                         ('Held', filters.Held),
                         ('Simoncelli', filters.Simoncelli),
                         ('Papadakis', filters.Papadakis)]:
    filt = FiltClass(g)
    bounds = filt.estimate_frame_bounds()
    kernel_evals = filt.evaluate(g.e)
    # Sum of squares should be ~1 for tight frames
    sos = sum(k**2 for k in kernel_evals)

    tight_data[name] = {
        'graph': 'Ring_16',
        'A': float(bounds[0]),
        'B': float(bounds[1]),
        'n_filters': len(kernel_evals),
        'kernel_evals': [to_list(k) for k in kernel_evals],
        'sum_of_squares': to_list(sos),
        'lmax': float(g.lmax),
    }

# Add tight frame data to frames file
frames_data['tight_frames'] = tight_data
save_golden('golden_frames.json', frames_data)


# =============================================================================
# M. Bank kernel evaluations: per-band kernel values at eigenvalues
#    for MexicanHat and Meyer (with rgsp-mapped equivalents)
# =============================================================================
print("M. Generating bank kernel evaluation golden values ...")

# M1: MexicanHat bank on Ring(20) -- per-band kernel evals already in filterbanks
# Add filtered output for a test signal
g = graphs.Ring(20)
g.compute_fourier_basis()
g.estimate_lmax()
Nf = 4
filt = filters.MexicanHat(g, Nf=Nf)

np.random.seed(42)
x_mh = np.random.randn(20)
y_mh = filt.filter(x_mh, method='exact')

# Store per-band filtered output
filterbank_data['mexicanhat_ring20']['signal'] = to_list(x_mh)
filterbank_data['mexicanhat_ring20']['filtered_output_per_band'] = [
    to_list(y_mh[:, i]) for i in range(y_mh.shape[1])
]

# M2: Meyer bank on Path(16) -- add filtered output
g = graphs.Path(16)
g.compute_fourier_basis()
g.estimate_lmax()
Nf = 3
filt = filters.Meyer(g, Nf=Nf)

np.random.seed(42)
x_my = np.random.randn(16)
y_my = filt.filter(x_my, method='exact')

filterbank_data['meyer_path16']['signal'] = to_list(x_my)
filterbank_data['meyer_path16']['filtered_output_per_band'] = [
    to_list(y_my[:, i]) for i in range(y_my.shape[1])
]

# M3: Heat bank at multiple scales on Ring(12)
g = graphs.Ring(12)
g.compute_fourier_basis()
g.estimate_lmax()
heat_scales = [0.1, 0.5, 2.0, 5.0]
heat_bank_kernels = []
for s in heat_scales:
    kernel_vals = np.exp(-s * g.e / g.lmax)
    heat_bank_kernels.append(to_list(kernel_vals))

np.random.seed(42)
x_hb = np.random.randn(12)
heat_bank_filtered = []
for s in heat_scales:
    filt = filters.Heat(g, scale=s)
    y = filt.filter(x_hb, method='exact')
    heat_bank_filtered.append(to_list(y.flatten()))

filterbank_data['heat_bank_ring12'] = {
    'graph': 'Ring_12',
    'scales': heat_scales,
    'lmax': float(g.lmax),
    'eigenvalues': to_list(np.sort(g.e)),
    'kernel_evals_per_scale': heat_bank_kernels,
    'signal': to_list(x_hb),
    'filtered_output_per_scale': heat_bank_filtered,
    'note': 'Heat kernel at multiple scales. PyGSP: exp(-scale*lambda/lmax), rgsp: exp(-t*lambda) with t=scale/lmax.',
}

save_golden('golden_filterbanks.json', filterbank_data)


# =============================================================================
# N. Classification Tikhonov: Ring(20) with 3-class labels
# =============================================================================
print("N. Generating classification_tikhonov golden values ...")

g = graphs.Ring(20)
g.compute_fourier_basis()
g.compute_laplacian()

# Create 3-class labels (deterministic, no randomness needed)
labels = np.zeros(20, dtype=int)
labels[0:7] = 0
labels[7:14] = 1
labels[14:20] = 2

tau = 5.0
mask = np.ones(g.N, dtype=bool)

# classification_tikhonov converts to logits and calls regression_tikhonov
# We replicate this manually to get the logits output
Y_logits = np.zeros((20, 3))
for i in range(3):
    Y_logits[:, i] = (labels == i).astype(float)

# Apply regression_tikhonov to each column
output_logits = np.zeros((20, 3))
for i in range(3):
    output_logits[:, i] = learning.regression_tikhonov(g, Y_logits[:, i], mask, tau)

# Predicted classes = argmax
predicted = np.argmax(output_logits, axis=1)

learning_data['classification_ring20'] = {
    'graph': 'Ring_20',
    'tau': tau,
    'labels': to_list(labels),
    'n_classes': 3,
    'output_logits': [to_list(output_logits[:, i]) for i in range(3)],
    'predicted_classes': to_list(predicted),
    'note': 'PyGSP classification_tikhonov: one-hot -> regression_tikhonov per column. rgsp: classification_tikhonov(g, labels, tau).',
}

save_golden('golden_learning.json', learning_data)


# =============================================================================
# O. Normalized Laplacian eigenvalues
# =============================================================================
print("O. Generating normalized Laplacian golden values ...")

norm_lap_data = {}
for name, constructor in [('Ring_10', lambda: graphs.Ring(10)),
                           ('Grid2d_3_4', lambda: graphs.Grid2d(3, 4))]:
    g = constructor()
    g.compute_fourier_basis()
    # Standard (combinatorial) eigenvalues already in golden_graphs.json
    # Now compute normalized Laplacian eigenvalues
    A = g.W
    d = np.array(A.sum(axis=1)).flatten()
    D_inv_sqrt = sparse.diags(1.0 / np.sqrt(np.maximum(d, 1e-10)))
    L_norm = sparse.eye(g.N) - D_inv_sqrt @ A @ D_inv_sqrt
    eigs_norm = np.sort(np.linalg.eigvalsh(L_norm.toarray()))

    norm_lap_data[name] = {
        'n': int(g.N),
        'eigenvalues_normalized': to_list(eigs_norm),
        'note': 'Normalized Laplacian: I - D^{-1/2} A D^{-1/2}. Eigenvalues in [0, 2].',
    }

save_golden('golden_normalized_laplacian.json', norm_lap_data)


# =============================================================================
# P. Lanczos filter: Ring(12) Heat(0.5) K=30 -- compare to exact
# =============================================================================
print("P. Generating Lanczos filter golden values ...")

g = graphs.Ring(12)
g.compute_fourier_basis()
g.estimate_lmax()
scale = 0.5
filt = filters.Heat(g, scale=scale)

np.random.seed(42)
x_lanc = np.random.randn(12)

# Exact filtered output (same as section C)
y_exact_lanc = filt.filter(x_lanc, method='exact')
# Chebyshev-filtered output at order 30
y_cheby_lanc = filt.filter(x_lanc, method='chebyshev', order=30)

lanczos_data = {
    'graph': 'Ring_12',
    'scale': scale,
    'lmax': float(g.lmax),
    'signal': to_list(x_lanc),
    'filtered_exact': to_list(y_exact_lanc.flatten()),
    'filtered_chebyshev_m30': to_list(y_cheby_lanc.flatten()),
    'note': 'Lanczos filtering should match exact output closely. rgsp filter_signal_lanczos(g, x, kernel, K=30).',
}

save_golden('golden_lanczos.json', lanczos_data)


# =============================================================================
# Q. SGWT round-trip: Ring(20) MexicanHat, forward then inverse
# =============================================================================
print("Q. Generating SGWT round-trip golden values ...")

g = graphs.Ring(20)
g.compute_fourier_basis()
g.estimate_lmax()

np.random.seed(42)
x_sgwt = np.random.randn(20)

# Forward transform with MexicanHat, Nf=3 (lowpass + 2 bandpass)
Nf = 3
filt_sgwt = filters.MexicanHat(g, Nf=Nf)
y_sgwt = filt_sgwt.filter(x_sgwt, method='exact')

# Per-band output
sgwt_bands = [to_list(y_sgwt[:, i]) for i in range(y_sgwt.shape[1])]

# Frame bounds for reconstruction quality assessment
bounds_sgwt = filt_sgwt.estimate_frame_bounds()

# Energy preservation: sum of band energies
band_energies = [float(np.sum(y_sgwt[:, i]**2)) for i in range(y_sgwt.shape[1])]
input_energy = float(np.sum(x_sgwt**2))

sgwt_data = {
    'graph': 'Ring_20',
    'Nf': Nf,
    'lmax': float(g.lmax),
    'eigenvalues': to_list(np.sort(g.e)),
    'signal': to_list(x_sgwt),
    'n_bands': y_sgwt.shape[1],
    'band_outputs': sgwt_bands,
    'input_energy': input_energy,
    'band_energies': band_energies,
    'frame_bounds': {'A': float(bounds_sgwt[0]), 'B': float(bounds_sgwt[1])},
    'scales': to_list(np.array(filt_sgwt.scales)),
    'note': 'SGWT forward: filter bank applied to signal. Energy should satisfy A*||x||^2 <= sum(||W_i x||^2) <= B*||x||^2.',
}

save_golden('golden_sgwt.json', sgwt_data)


print("\nDone! All golden values generated.")
