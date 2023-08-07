# Requirements

The API for this library needs to be flexible enough to be used for a variety of applications.
As a handy guide, the following paper outlines some requirements:
- Barrett, John W., and L Glaser. ‘Monte Carlo Simulations of Random Non-Commutative Geometries’. Journal of Physics A: Mathematical and Theoretical 49, no. 24 (17 June 2016): 245001. https://doi.org/10.1088/1751-8113/49/24/245001.

They are summarised below:
1. The actions need to be 'spectral': we need $S(D) = \sum_i V(\lambda_i)$ for some potential function $V$, such that $V\geq b$ for some $b\in \mathbb{R}$, and $\lambda_i$ are the eigenvalues of the Dirac operator $D$.
1. The actions need to asymptote to infinite: $S(D) = \sum V(D)$, then $V(x) \to \infty$ as $x \to \infty$.
1. The ability to make measurements of 'observables'. Observables are functions of the Dirac operator. For Monte Carlo simulations, this boils down to formulas of this: $\langle f \rangle_N = \frac{1}{N} \sum_{j=1}^N f(D_j)$, where $\{D_j\}$ is an ensemble of Dirac operators sampled from the simulation.
1. The integrated autocorrelation of measured observables needs to be calculated.
