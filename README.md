# Porteous Meets Walsh

Code accompanying

> S. Lagauw, L. Kayser, and B. De Moor, *Porteous Meets Walsh: Counting Critical Points for $\mathcal{H}_2$-Norm Model Order Reduction Using Algebraic Geometry Behind Multiparameter Eigenvalue Problems*, 2026. [arXiv:2611.00258](https://arxiv.org/abs/2611.00258)

## Abstract

Several polynomial systems arising in systems and control have been shown to exhibit a partially linear structure, which allows them to be reformulated as a multiparameter eigenvalue problem (MEP). In this work, we show how such MEP formulations can be leveraged to obtain refined information about the solution set of the polynomial system. We demonstrate, in the setting of single-input single-output $\mathcal{H}_2$-norm model order reduction, how the Thom–Porteous formula can be used to derive an exact count of the complex solutions of the critical equations derived from Walsh's theorem. We prove, under suitable genericity assumptions, that their set of complex solutions is finite and that all solutions occur without multiplicity. Moreover, we derive an explicit formula for the total number of complex solutions as a function of the model orders. Throughout, we provide a self-contained and application-oriented exposition of the algebro-geometric tools required for our analysis, with the aim of making them accessible to the broader engineering community.

For a generic model of McMillan degree $N$, reduction to order $1 \le n < N$ has $\sum_{j=0}^{n} \binom{N-n-1+j}{j}\, 2^j$ critical points (Theorem 3.10).

## Repository structure

```
MATLAB/
  - getMep.m, getRestrictedMep.m      MEP formulations of the critical equations
  - getOptimalZeros.m, h2Error.m      optimal numerator and H2 error for a fixed denominator
  - expectedDegree.m                  count of Theorem 3.10
  - symToMep.m, tfCoefficients.m,
    getToeplitz.m                     internal helpers
  - example_5_*.m                     examples of Section 5
m2/
  - WalshVariety.m2                   Walsh variety, MEP and bad locus as ideals
  - Examples.m2                       Examples 4.1, 5.1, 5.3, 5.4, 5.5 and 5.8
```

In MATLAB, MEPs are solved numerically with MacaulayLab, whereas Macaulay2 code computes dimensions, degrees and bad loci from the defining ideals.

## Examples

| Paper | Files | Content |
|---|---|---|
| Ex. 4.1 | `Examples.m2` | Bad loci for denominators with double roots. |
| Ex. 5.1 | `example_5_1.m`, `Examples.m2` | Continuous time, $N=3$, $n=1$: five critical points, two local minimizers. |
| Ex. 5.2 | `example_5_2.m` | Discrete time, $N=6$, $n=2$: 49 critical points; contour plot of Figure 1. |
| Ex. 5.3 | `example_5_3.m`, `Examples.m2` | Double pole, $n=1$: a bad locus of one point leaves 4 critical points. |
| Ex. 5.4 | `example_5_4.m`, `Examples.m2` | Same model, $n=2$: a double bad-locus point leaves 5 critical points. |
| Ex. 5.5 | `example_5_5.m`, `Examples.m2` | $N=2$, $n=1$: a critical point of multiplicity two. |
| Ex. 5.6 | `example_5_6.m` | Model of Ex. 5.2 restricted to a double pole: polynomial eigenvalue problem $M'(p)$. |
| Ex. 5.8 | `example_5_8.m`, `Examples.m2` | $N=4$, $n=3$: 15 critical points; the global minimizer $1/(s+1)^3$ has a triple pole. |

## Functions

**MATLAB.** The model is a SISO `tf` object `sys`, continuous-time if `sys.Ts == 0`. Polynomials are row vectors in descending powers.

- `[mep, eqs, M] = getMep(sys, n)`: the MEP (9) for reduction to order `n`, with the critical equations (6) and the matrix $M(a,c,d)$. Each eigentuple $(a_1,\dots,a_n)$ gives the denominator `[1, a1, ..., an]`; the bad locus is described by `M(:, 2:end)`.
- `getRestrictedMep(sys, n)`: the same, restricted to a single pole of multiplicity `n` (Examples 5.6 and 5.8).
- `b = getOptimalZeros(sys, a)` and `J = h2Error(sys, a)`: the optimal numerator (21) and the $\mathcal{H}_2$ error for a fixed stable denominator `a`.
- `expectedDegree(N, n)`: the count (14).

**Macaulay2** (`WalshVariety.m2`; run `Examples.m2` from inside `m2/` so that `load "WalshVariety.m2"` finds the file). The functions take the orders `(N, n)`, the option `Discrete` (default `false`) and specific data `CD => (c, d)` (`C => c` for `badLocus`).

- `Walsh(N, n)`: ideal of the Walsh variety $\mathcal{W}$ (Definition 3.1), in the chart $a_n = c_N = 1$.
- `Mmatrix(N, n)` and `MEP(N, n)`: the matrix $M(a,c,d)$ and the ideal of its maximal minors, which defines $\mathcal{D}$.
- `badLocus(N, n)`: ideal of the bad locus $\mathcal{B}$.
- `expectedDegree(N, n)`: the count (14).

## Requirements

- MATLAB with the Control System and Symbolic Math Toolboxes, and [MacaulayLab](https://www.macaulaylab.net).
- [Macaulay2](https://macaulay2.com) with the packages `NumericSolutions` and `RealRoots`.

A Julia implementation of the MEP formulation is available at <https://codeocean.com/capsule/4036835>.

## Citation

```bibtex
@misc{LagauwKayserDeMoor2026,
  author        = {Lagauw, Sibren and Kayser, Leonie and De Moor, Bart},
  title         = {Porteous Meets {Walsh}: Counting Critical Points for $\mathcal{H}_2$-Norm Model Order Reduction Using Algebraic Geometry Behind Multiparameter Eigenvalue Problems},
  year          = {2026},
  eprint        = {2611.00258},
  archivePrefix = {arXiv},
  url           = {https://arxiv.org/abs/2611.00258}
}
```
