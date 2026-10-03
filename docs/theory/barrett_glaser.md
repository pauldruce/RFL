# The Barrett-Glaser Spectral Action

This document provides the mathematical derivations, trace factorisations, and local variation formulae implemented by `rfl::BarrettGlaserAction` in `src/core/BarrettGlaser/BarrettGlaserAction.hpp`.

---

## 1. Noncommutative Geometry & Finite Spectral Triples

A finite spectral triple $(\mathcal{A}, \mathcal{H}, D, \gamma, J)$ models a fuzzy space in noncommutative geometry.

In the Barrett-Glaser formulation, the Hilbert space is the tensor product of a spinor space and the algebra of complex matrices:

$$
\mathcal{H} = \mathbb{C}^{d_\gamma} \otimes \mathcal{M}_N(\mathbb{C}) \cong \mathbb{C}^{d_\gamma} \otimes \mathbb{C}^N \otimes \mathbb{C}^N
$$

where:
* $d_\gamma$ is the dimension of the Clifford algebra representation.
* $N$ is the matrix dimension.
* $V = p + q$ is the number of Clifford generators.

The Dirac operator $D: \mathcal{H} \to \mathcal{H}$ acts on a state $\psi \in \mathcal{H}$ as:

$$
D(\psi) = \sum_{a=1}^V \gamma^a \otimes (M_a \psi + \epsilon_a \psi M_a)
$$

In Kronecker product notation on $\mathbb{C}^{d_\gamma} \otimes \mathbb{C}^N \otimes \mathbb{C}^N$:

$$
D = \sum_{a=1}^V \gamma^a \otimes \mathcal{D}_a, \quad \text{with } \mathcal{D}_a = M_a \otimes \mathbf{1}_N + \epsilon_a \mathbf{1}_N \otimes M_a^T
$$

where $\epsilon_a \in \{ +1, -1 \}$ are the Clifford generator signs: $\epsilon_a = +1$ for $a \le p$ (Hermitian) and $\epsilon_a = -1$ for $a > p$ (anti-Hermitian).

The gamma matrices satisfy the Clifford anticommutation relations:

$$
\{ \gamma^a, \gamma^b \} = 2 \epsilon_a \delta^{ab} \mathbf{1}_{d_\gamma}
$$

and their spinor trace satisfies orthogonality:

$$
\mathrm{tr}(\gamma^a \gamma^b) = d_\gamma \epsilon_a \delta^{ab}
$$

---

## 2. The Spectral Action Functional

The Barrett-Glaser action is defined by a polynomial functional of the Dirac operator:

$$
S(D) = g_2 \mathrm{Tr}_{\mathcal{H}}(D^2) + g_4 \mathrm{Tr}_{\mathcal{H}}(D^4)
$$

where $g_2$ and $g_4$ are coupling constants.

### 2.1 Quadratic Trace $\mathrm{Tr}(D^2)$

Expanding $D^2$ and applying Clifford trace orthogonality $\mathrm{tr}(\gamma^a \gamma^b) = d_\gamma \epsilon_a \delta^{ab}$:

$$
\mathrm{Tr}_{\mathcal{H}}(D^2) = \sum_{a,b=1}^V \mathrm{tr}(\gamma^a \gamma^b) \cdot \mathrm{Tr}_{\mathrm{alg}}(\mathcal{D}_a \mathcal{D}_b) = d_\gamma \sum_{a=1}^V \epsilon_a \mathrm{Tr}_{\mathrm{alg}}(\mathcal{D}_a^2)
$$

Expanding the algebra operator $\mathcal{D}_a = M_a \otimes \mathbf{1}_N + \epsilon_a \mathbf{1}_N \otimes M_a^T$:

$$
\mathcal{D}_a^2 = M_a^2 \otimes \mathbf{1}_N + 2 \epsilon_a M_a \otimes M_a^T + \mathbf{1}_N \otimes (M_a^T)^2
$$

Taking the trace over $\mathbb{C}^N \otimes \mathbb{C}^N$ yields:

$$
\mathrm{Tr}_{\mathrm{alg}}(\mathcal{D}_a^2) = 2 N \mathrm{Tr}(M_a^2) + 2 \epsilon_a (\mathrm{Tr}(M_a))^2
$$

Multiplying by $d_\gamma \epsilon_a$ gives the component formula:

$$
\mathrm{Tr}_{\mathcal{H}}(D^2) = 2 d_\gamma \sum_{a=1}^V \left[ \epsilon_a N \mathrm{Tr}(M_a^2) + (\mathrm{Tr}(M_a))^2 \right]
$$

### 2.2 Quartic Trace $\mathrm{Tr}(D^4)$

Expanding $D^4$ across the tensor product space:

$$
\mathrm{Tr}_{\mathcal{H}}(D^4) = \sum_{a,b,c,d=1}^V \Omega_{abcd} \cdot \mathrm{Tr}_{\mathrm{alg}}(\mathcal{D}_a \mathcal{D}_b \mathcal{D}_c \mathcal{D}_d)
$$

where $\Omega_{abcd} = \mathrm{tr}(\gamma^a \gamma^b \gamma^c \gamma^d)$ is the 4-gamma Clifford trace tensor.

The contractions partition into three symmetry classes:
1. **$B_4$:** Terms where all four indices $(a, b, c, d)$ are distinct.
2. **$B_2$:** Terms with two distinct index pairs (e.g. $\Omega_{aabb}, \Omega_{abab}, \Omega_{abba}$).
3. **$B$:** Quartic self-couplings where $a = b = c = d$.

---

## 3. Elementary Perturbations & Local Variations

In MCMC simulations, a local move perturbs a single matrix element of component matrix $M_x$:

$$
M_x \to M_x + \delta M
$$

To maintain Hermiticity of $M_x$, an update at row $i$ and column $j$ with complex parameter $z$ takes the form:

$$
\delta M = z E_{ij} + \bar{z} E_{ji}
$$

where $E_{ij}$ is the standard matrix unit: $(E_{ij})_{kl} = \delta_{ik} \delta_{jl}$.

### 3.1 Diagonal vs Off-Diagonal Perturbations

**Off-diagonal move ($i \neq j$):**
Updates both $(i, j)$ and $(j, i)$:

$$
\delta M_{ij} = z, \quad \delta M_{ji} = \bar{z}
$$

The trace vanishes: $\mathrm{Tr}(\delta M) = 0$.

**Diagonal move ($i = j$):**
The variation collapses to:

$$
\delta M_{ii} = z + \bar{z} = 2 \mathrm{Re}(z)
$$

The imaginary component $\mathrm{Im}(z)$ has zero effect, ensuring $\mathrm{Im}(M_{ii}) = 0$.

---

### 3.2 Quadratic Variation $\Delta S_2$: From $D$ to $M$

The exact change in the quadratic trace is:

$$
\Delta \mathrm{Tr}(D^2) = \mathrm{Tr}_{\mathcal{H}}((D + \delta D)^2) - \mathrm{Tr}_{\mathcal{H}}(D^2)
$$

**Step 1 (Operator expansion):**
Expand the square and apply trace cyclicity $\mathrm{Tr}_{\mathcal{H}}(D \delta D) = \mathrm{Tr}_{\mathcal{H}}(\delta D D)$:

$$
\Delta \mathrm{Tr}(D^2) = 2 \mathrm{Tr}_{\mathcal{H}}(D \delta D) + \mathrm{Tr}_{\mathcal{H}}((\delta D)^2)
$$

**Step 2 (Clifford trace factorisation):**
Substitute $D = \sum_{a=1}^V \gamma^a \otimes \mathcal{D}_a$ and $\delta D = \gamma^x \otimes \delta \mathcal{D}_x$.
Factor the Hilbert space trace $\mathrm{Tr}_{\mathcal{H}} = \mathrm{tr} \otimes \mathrm{Tr}_{\mathrm{alg}}$:

$$
\mathrm{Tr}_{\mathcal{H}}(D \delta D) = \sum_{a=1}^V \mathrm{tr}(\gamma^a \gamma^x) \cdot \mathrm{Tr}_{\mathrm{alg}}(\mathcal{D}_a \delta \mathcal{D}_x)
$$

Apply Clifford trace orthogonality $\mathrm{tr}(\gamma^a \gamma^x) = d_\gamma \epsilon_x \delta^{ax}$. Only the $x$-th term survives:

$$
\mathrm{Tr}_{\mathcal{H}}(D \delta D) = d_\gamma \epsilon_x \mathrm{Tr}_{\mathrm{alg}}(\mathcal{D}_x \delta \mathcal{D}_x)
$$

Similarly, for the variation squared:

$$
\mathrm{Tr}_{\mathcal{H}}((\delta D)^2) = \mathrm{tr}((\gamma^x)^2) \cdot \mathrm{Tr}_{\mathrm{alg}}((\delta \mathcal{D}_x)^2) = d_\gamma \epsilon_x \mathrm{Tr}_{\mathrm{alg}}((\delta \mathcal{D}_x)^2)
$$

Combining these yields the algebra-level variation:

$$
\Delta \mathrm{Tr}(D^2) = d_\gamma \epsilon_x \left[ 2 \mathrm{Tr}_{\mathrm{alg}}(\mathcal{D}_x \delta \mathcal{D}_x) + \mathrm{Tr}_{\mathrm{alg}}((\delta \mathcal{D}_x)^2) \right]
$$

**Step 3 (Expansion into matrix components $M_x$ and $\delta M$):**
Substitute $\mathcal{D}_x = M_x \otimes \mathbf{1}_N + \epsilon_x \mathbf{1}_N \otimes M_x^T$ and $\delta \mathcal{D}_x = \delta M \otimes \mathbf{1}_N + \epsilon_x \mathbf{1}_N \otimes \delta M^T$.
Expand the operator product on $\mathbb{C}^N \otimes \mathbb{C}^N$:

$$
\mathcal{D}_x \delta \mathcal{D}_x = (M_x \delta M) \otimes \mathbf{1}_N + \epsilon_x M_x \otimes \delta M^T + \epsilon_x \delta M \otimes M_x^T + \mathbf{1}_N \otimes (M_x \delta M)^T
$$

Evaluate the trace using $\mathrm{Tr}(\mathbf{1}_N) = N$ and $\mathrm{Tr}(A^T) = \mathrm{Tr}(A)$:

$$
\mathrm{Tr}_{\mathrm{alg}}(\mathcal{D}_x \delta \mathcal{D}_x) = 2 N \mathrm{Tr}(M_x \delta M) + 2 \epsilon_x \mathrm{Tr}(M_x) \mathrm{Tr}(\delta M)
$$

Setting $M_x \to \delta M$ gives the square term:

$$
\mathrm{Tr}_{\mathrm{alg}}((\delta \mathcal{D}_x)^2) = 2 N \mathrm{Tr}((\delta M)^2) + 2 \epsilon_x (\mathrm{Tr}(\delta M))^2
$$

Multiply by $d_\gamma \epsilon_x$ (using $\epsilon_x^2 = 1$) to obtain the component matrix form:

$$
\Delta \mathrm{Tr}(D^2) = 4 d_\gamma \left[ \epsilon_x N \mathrm{Tr}(M_x \delta M) + \mathrm{Tr}(M_x) \mathrm{Tr}(\delta M) \right] + 2 d_\gamma \left[ \epsilon_x N \mathrm{Tr}((\delta M)^2) + (\mathrm{Tr}(\delta M))^2 \right]
$$

**Step 4 (Evaluation on elementary move $\delta M = z E_{ij} + \bar{z} E_{ji}$):**
Use $\mathrm{Tr}(M_x \delta M) = 2 \mathrm{Re}(z M_x(j, i))$:

**Off-diagonal move ($i \neq j$):**
Here $\mathrm{Tr}(\delta M) = 0$ and $\mathrm{Tr}((\delta M)^2) = 2 |z|^2$:

$$
\Delta_2 = 4 d_\gamma N \left( 2 \mathrm{Re}(z M_x(j, i)) + |z|^2 \right)
$$

**Diagonal move ($i = j$):**
Here $\mathrm{Tr}(\delta M) = 2 \mathrm{Re}(z)$ and $\mathrm{Tr}((\delta M)^2) = 4 \mathrm{Re}(z)^2$:

$$
\Delta_2 = 8 d_\gamma \mathrm{Re}(z) \left[ N (M_x(i, i) + \mathrm{Re}(z)) + \epsilon_x (\mathrm{Tr}(M_x) + \mathrm{Re}(z)) \right]
$$

---

### 3.3 Quartic Variation $\Delta S_4$: From $D$ to $M$

The exact change in the quartic trace is:

$$
\Delta \mathrm{Tr}(D^4) = \mathrm{Tr}_{\mathcal{H}}((D + \delta D)^4) - \mathrm{Tr}_{\mathcal{H}}(D^4)
$$

**Step 1 (Operator binomial expansion):**
Expand $((D + \delta D)^4 - D^4)$ and apply trace cyclicity to group terms by perturbation order:

$$
\Delta \mathrm{Tr}(D^4) = \underbrace{4 \mathrm{Tr}_{\mathcal{H}}(D^3 \delta D)}_{\text{Linear } O(z)} + \underbrace{2 \mathrm{Tr}_{\mathcal{H}}(D^2 (\delta D)^2) + \mathrm{Tr}_{\mathcal{H}}(D \delta D D \delta D)}_{\text{Quadratic } O(z^2)} + \underbrace{4 \mathrm{Tr}_{\mathcal{H}}(D (\delta D)^3)}_{\text{Cubic } O(z^3)} + \underbrace{\mathrm{Tr}_{\mathcal{H}}((\delta D)^4)}_{\text{Quartic } O(z^4)}
$$

**Step 2 (Clifford trace factorisation):**
Substitute $D = \sum_{a=1}^V \gamma^a \otimes \mathcal{D}_a$ and $\delta D = \gamma^x \otimes \delta \mathcal{D}_x$.
Factoring $\mathrm{Tr}_{\mathcal{H}} = \mathrm{tr} \otimes \mathrm{Tr}_{\mathrm{alg}}$ isolates the Clifford trace tensors $\Omega_{abcd} = \mathrm{tr}(\gamma^a \gamma^b \gamma^c \gamma^d)$:

**Linear term ($O(z)$):**

$$
4 \mathrm{Tr}_{\mathcal{H}}(D^3 \delta D) = 4 \sum_{i_1, i_2, i_3=1}^V \Omega_{i_1 i_2 i_3 x} \cdot \mathrm{Tr}_{\mathrm{alg}}(\mathcal{D}_{i_1} \mathcal{D}_{i_2} \mathcal{D}_{i_3} \delta \mathcal{D}_x)
$$

**Quadratic terms ($O(z^2)$):**

$$
\begin{aligned}
2 \mathrm{Tr}_{\mathcal{H}}(D^2 (\delta D)^2) &= 2 d_\gamma \epsilon_x \sum_{a=1}^V \epsilon_a \mathrm{Tr}_{\mathrm{alg}}(\mathcal{D}_a^2 (\delta \mathcal{D}_x)^2) \\
\mathrm{Tr}_{\mathcal{H}}(D \delta D D \delta D) &= \sum_{a,b=1}^V \Omega_{axbx} \cdot \mathrm{Tr}_{\mathrm{alg}}(\mathcal{D}_a \delta \mathcal{D}_x \mathcal{D}_b \delta \mathcal{D}_x)
\end{aligned}
$$

**Cubic term ($O(z^3)$):**
Since $(\gamma^x)^3 = \epsilon_x \gamma^x$, Clifford orthogonality leaves only $a = x$:

$$
4 \mathrm{Tr}_{\mathcal{H}}(D (\delta D)^3) = 4 d_\gamma \mathrm{Tr}_{\mathrm{alg}}(\mathcal{D}_x (\delta \mathcal{D}_x)^3)
$$

**Quartic term ($O(z^4)$):**
Since $(\gamma^x)^4 = \mathbf{1}_{d_\gamma}$:

$$
\mathrm{Tr}_{\mathcal{H}}((\delta D)^4) = d_\gamma \mathrm{Tr}_{\mathrm{alg}}((\delta \mathcal{D}_x)^4)
$$

**Step 3 (Expansion into matrix components and algebra traces):**

Each algebra operator has the Kronecker form $\mathcal{D}_a = M_a \otimes \mathbf{1}_N + \epsilon_a \mathbf{1}_N \otimes M_a^T$.
Expanding operator products on $\mathbb{C}^N \otimes \mathbb{C}^N$ and taking the trace yields the component matrix expressions.

#### 1. Linear variation $\mathrm{Tr}_{\mathrm{alg}}(\mathcal{D}_{i_1} \mathcal{D}_{i_2} \mathcal{D}_{i_3} \delta \mathcal{D}_x)$

Substituting $\mathcal{D}_a = M_a \otimes \mathbf{1}_N + \epsilon_a \mathbf{1}_N \otimes M_a^T$ into the product and taking the algebra trace $\mathrm{Tr}_{\mathrm{alg}} = \mathrm{Tr} \otimes \mathrm{Tr}$ contracts the 16 tensor products into three trace topologies:

**Single-trace (leading $O(N)$):**

$$
N \left[ \mathrm{Tr}(M_1 M_2 M_3 \delta M) + \epsilon_1 \epsilon_2 \epsilon_3 \epsilon_x \mathrm{Tr}(M_3 M_2 M_1 \delta M) \right]
$$

**Double-trace (matrix-matrix and matrix-perturbation):**

$$
\begin{aligned}
& \left[ \epsilon_3 \mathrm{Tr}(M_1 M_2 \delta M) + \epsilon_1 \epsilon_2 \epsilon_x \mathrm{Tr}(M_2 M_1 \delta M) \right] \mathrm{Tr}(M_3) \\
+\;& \left[ \epsilon_2 \mathrm{Tr}(M_1 M_3 \delta M) + \epsilon_1 \epsilon_3 \epsilon_x \mathrm{Tr}(M_3 M_1 \delta M) \right] \mathrm{Tr}(M_2) \\
+\;& \left[ \epsilon_1 \mathrm{Tr}(M_2 M_3 \delta M) + \epsilon_2 \epsilon_3 \epsilon_x \mathrm{Tr}(M_3 M_2 \delta M) \right] \mathrm{Tr}(M_1) \\
+\;& (\epsilon_1 \epsilon_2 + \epsilon_3 \epsilon_x) \mathrm{Tr}(M_1 M_2) \mathrm{Tr}(M_3 \delta M) \\
+\;& (\epsilon_1 \epsilon_3 + \epsilon_2 \epsilon_x) \mathrm{Tr}(M_1 M_3) \mathrm{Tr}(M_2 \delta M) \\
+\;& (\epsilon_2 \epsilon_3 + \epsilon_1 \epsilon_x) \mathrm{Tr}(M_2 M_3) \mathrm{Tr}(M_1 \delta M)
\end{aligned}
$$

**Trace-perturbation:**

$$
\left[ \epsilon_x \mathrm{Tr}(M_1 M_2 M_3) + \epsilon_1 \epsilon_2 \epsilon_3 \mathrm{Tr}(M_3 M_2 M_1) \right] \mathrm{Tr}(\delta M)
$$

These terms map directly to variables `t_1` to `t_8` in `BarrettGlaserAction::delta4`.

#### 2. Quadratic variations $\mathrm{Tr}_{\mathrm{alg}}(\mathcal{D}_a^2 (\delta \mathcal{D}_x)^2)$ and $\mathrm{Tr}_{\mathrm{alg}}(\mathcal{D}_a \delta \mathcal{D}_x \mathcal{D}_a \delta \mathcal{D}_x)$

Expand the operator squares:

$$
\begin{aligned}
\mathcal{D}_a^2 &= M_a^2 \otimes \mathbf{1}_N + 2 \epsilon_a M_a \otimes M_a^T + \mathbf{1}_N \otimes (M_a^T)^2 \\
(\delta \mathcal{D}_x)^2 &= (\delta M)^2 \otimes \mathbf{1}_N + 2 \epsilon_x \delta M \otimes \delta M^T + \mathbf{1}_N \otimes (\delta M^T)^2
\end{aligned}
$$

Evaluate the algebra traces:

$$
\begin{aligned}
\mathrm{Tr}_{\mathrm{alg}}(\mathcal{D}_a^2 (\delta \mathcal{D}_x)^2) &= 2 N \mathrm{Tr}(M_a^2 (\delta M)^2) + 4 \epsilon_a \epsilon_x (\mathrm{Tr}(M_a \delta M))^2 + 2 \mathrm{Tr}(M_a^2) \mathrm{Tr}((\delta M)^2) \\
\mathrm{Tr}_{\mathrm{alg}}(\mathcal{D}_a \delta \mathcal{D}_x \mathcal{D}_a \delta \mathcal{D}_x) &= 2 N \mathrm{Tr}((M_a \delta M)^2) + 4 \epsilon_a \epsilon_x \mathrm{Tr}(M_a^2) \mathrm{Tr}((\delta M)^2) + 2 (\mathrm{Tr}(M_a \delta M))^2
\end{aligned}
$$

#### 3. Cubic variation $\mathrm{Tr}_{\mathrm{alg}}(\mathcal{D}_x (\delta \mathcal{D}_x)^3)$

Expand $(\delta \mathcal{D}_x)^3$:

$$
(\delta \mathcal{D}_x)^3 = (\delta M)^3 \otimes \mathbf{1}_N + 3 \epsilon_x (\delta M)^2 \otimes \delta M^T + 3 \delta M \otimes (\delta M^T)^2 + \epsilon_x \mathbf{1}_N \otimes (\delta M^T)^3
$$

Multiplying by $\mathcal{D}_x = M_x \otimes \mathbf{1}_N + \epsilon_x \mathbf{1}_N \otimes M_x^T$ and taking the trace yields:

$$
\mathrm{Tr}_{\mathrm{alg}}(\mathcal{D}_x (\delta \mathcal{D}_x)^3) = 2 N \mathrm{Tr}(M_x (\delta M)^3) + 6 \mathrm{Tr}(M_x \delta M) \mathrm{Tr}((\delta M)^2) + 2 \epsilon_x \mathrm{Tr}(M_x) \mathrm{Tr}((\delta M)^3)
$$

#### 4. Quartic variation $\mathrm{Tr}_{\mathrm{alg}}((\delta \mathcal{D}_x)^4)$

Expand the fourth power using the binomial theorem:

$$
(\delta \mathcal{D}_x)^4 = (\delta M)^4 \otimes \mathbf{1}_N + 4 \epsilon_x (\delta M)^3 \otimes \delta M^T + 6 (\delta M)^2 \otimes (\delta M^T)^2 + 4 \epsilon_x \delta M \otimes (\delta M^T)^3 + \mathbf{1}_N \otimes (\delta M^T)^4
$$

Evaluate the trace over $\mathbb{C}^N \otimes \mathbb{C}^N$:

$$
\mathrm{Tr}_{\mathrm{alg}}((\delta \mathcal{D}_x)^4) = 2 N \mathrm{Tr}((\delta M)^4) + 8 \epsilon_x \mathrm{Tr}((\delta M)^3) \mathrm{Tr}(\delta M) + 6 (\mathrm{Tr}((\delta M)^2))^2
$$

**Step 4 (Evaluation on elementary move $\delta M = z E_{ij} + \bar{z} E_{ji}$ in $O(N^2)$):**

Substituting the rank-1 matrix unit perturbation avoids full matrix multiplications:

Traces against $\delta M$ evaluate at indices $(i, j)$:

$$
\mathrm{Tr}(P \delta M) = z P(j, i) + \bar{z} P(i, j)
$$

where $P$ denotes precomputed matrix products ($M_{i_1} M_{i_2} M_{i_3}$, $M_{i_1} M_{i_2}$, etc.).

Higher powers contract algebraically via $(\delta M)^2 = |z|^2 (E_{ii} + E_{jj}) + (z^2 E_{ij}^2 + \bar{z}^2 E_{ji}^2)$.

The scalar quartic term evaluates to:

$$
\mathrm{Tr}_{\mathcal{H}}((\delta D)^4) = 4 d_\gamma (N + \epsilon_x) |z|^4
$$

Precomputing matrix products and evaluating trace updates at index $(i, j)$ reduces computational complexity from $O(d_\gamma^3 N^6)$ to $O(N^2)$ operations.

### 3.4 Total Variation & Acceptance

The total action variation is:

$$
\Delta S = g_2 \Delta_2 + g_4 \Delta_4
$$

The Metropolis acceptance probability is:

$$
P = \min\left(1, e^{-\Delta S}\right)
$$

---

## 4. Analytical Gradients & HMC Forces

In Hybrid Monte Carlo (HMC) simulations, conjugate momenta $P_k$ evolve along Hamilton's equations:

$$
\frac{d P_k}{d t} = -\nabla_{M_k} S(D) = -\left( g_2 \nabla_{M_k} \mathrm{Tr}(D^2) + g_4 \nabla_{M_k} \mathrm{Tr}(D^4) \right)
$$

### 4.1 Quadratic Gradient

Differentiating $\Delta_2$ with respect to matrix $M_k$ yields:

$$
\nabla_{M_k} \mathrm{Tr}(D^2) = \frac{\partial}{\partial M_k} \mathrm{Tr}(D^2) = 4 d_\gamma \left( N M_k + \epsilon_k \mathrm{Tr}(M_k) \mathbf{1}_N \right)
$$

### 4.2 Quartic Gradient

The quartic gradient partitions across the three index topologies:

$$
\nabla_{M_k} \mathrm{Tr}(D^4) = B_4(k) + B_2(k) + B(k)
$$

To preserve the Hermitian tangent space along numerical integrator trajectories, the force is projected onto the Hermitian subspace:

$$
\nabla_{\mathrm{Herm}} = 2 (\nabla + \nabla^\dagger)
$$

The factor of 2 accounts for the symmetric variation $\delta M = z E_{ij} + \bar{z} E_{ji}$.

---

## 5. References

* **[Barrett2016]:** Barrett, J. W. & Glaser, L. (2016). *Monte Carlo simulations of random non-commutative geometries.* Journal of Physics A: Mathematical and Theoretical, 49(24), 245001. [DOI: 10.1088/1751-8113/49/24/245001](https://doi.org/10.1088/1751-8113/49/24/245001) | [arXiv:1510.01377](https://arxiv.org/abs/1510.01377).
