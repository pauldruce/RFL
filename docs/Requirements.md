# Scientific Requirements

The library API must provide flexibility for diverse research applications. The following paper outlines core mathematical requirements:
* Barrett, John W., and Lisa Glaser. "Monte Carlo Simulations of Random Non-Commutative Geometries." *Journal of Physics A: Mathematical and Theoretical* 49, no. 24 (2016): 245001. [https://doi.org/10.1088/1751-8113/49/24/245001](https://doi.org/10.1088/1751-8113/49/24/245001).

Core requirements include:

1. **Spectral Actions:** Actions must be spectral: $S(D) = \sum_i V(\lambda_i)$ with bounded potential $V \ge b$ ($b \in \mathbb{R}$) and Dirac eigenvalues $\lambda_i$.
2. **Asymptotic Growth:** The potential function must grow asymptotically to infinity: $V(x) \to \infty$ as $x \to \infty$.
3. **Observable Measurement:** The library must measure observables $f(D)$. Monte Carlo simulations evaluate ensemble averages $\langle f \rangle_N = \frac{1}{N} \sum_{j=1}^N f(D_j)$ over sampled Dirac operators $\{D_j\}$.
4. **Autocorrelation Analysis:** The library must calculate the integrated autocorrelation time for measured observables.
