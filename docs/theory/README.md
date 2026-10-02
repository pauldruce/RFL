# Mathematical Theory & Derivations

This directory contains the mathematical derivations, physical invariants, and algorithmic factorisations implemented across the Random Fuzzy Library (RFL).

It serves as the authoritative Single Source of Truth (SSOT) bridging academic literature, the Obsidian Research Vault, and the C++20 / Python codebase.

---

## 1. Code-to-Theory Map

| Theory Note | Primary Scope & Derivations | Implementing C++ Classes | Python Bindings |
| :--- | :--- | :--- | :--- |
| **[Barrett-Glaser Action](barrett_glaser.md)** | Spectral action functional, local trace variations ($\Delta S_2, \Delta S_4$), and HMC force gradients. | `BarrettGlaserAction`, `Metropolis`, `Hamiltonian` | `rfl.Action`, `rfl.Metropolis` |

---

## 2. Global Mathematical Notation

All theory notes and code comments adhere to these mathematical symbols:

| Symbol | Concept | Description |
| :--- | :--- | :--- |
| $N$ | Matrix dimension | Dimension of the Hermitian matrices $M_a \in \mathcal{M}_N(\mathbb{C})$. |
| $(p, q)$ | Clifford signature / Type | The signature of the Clifford algebra $\mathcal{C}\ell(p, q)$, defining the spectral triple type $(p, q)$ with $p$ Hermitian and $q$ anti-Hermitian matrix generators. |
| $V$ | Generator count | Total number of Clifford generators: $V = p + q$. |
| $d_\gamma$ | Spinor dimension | Representation dimension of the Clifford algebra $\mathcal{C}\ell(p, q)$. |
| $\gamma^a$ | Gamma matrices | Generators satisfying $\{ \gamma^a, \gamma^b \} = 2 \epsilon_a \delta^{ab} \mathbf{1}_{d_\gamma}$. |
| $\epsilon_a$ | Clifford signs | Generator signs: $(\gamma^a)^2 = \epsilon_a \mathbf{1}_{d_\gamma}$, with $\epsilon_a = +1$ for $a \le p$ and $\epsilon_a = -1$ for $a > p$. |
| $M_a$ | Component matrices | Set of $V$ Hermitian matrices parameterising the geometry. |
| $D$ | Dirac operator | Assembled operator: $D = \sum_{a=1}^V \gamma^a \otimes M_a \in \mathcal{M}_{d_\gamma N}(\mathbb{C})$. |
| $g_2, g_4$ | Coupling constants | Quadratic ($g_2$) and quartic ($g_4$) coupling parameters. |
| $S(D)$ | Spectral action | Action functional: $S(D) = g_2 \mathrm{Tr}(D^2) + g_4 \mathrm{Tr}(D^4)$. |

---

## 3. Theory Note Blueprint

Every theory document in this directory follows this four-section structure:

1. **Mathematical Formulation:** Continuous and algebraic definitions, Hilbert spaces, and continuous symmetries.
2. **Discrete & Computational Realisation:** Tensor product factorisations, matrix algorithms, and computational complexity ($O(N^2)$ vs $O(N^3)$).
3. **Invariants & Symmetries:** Physical conservation laws and symmetry preservation (Hermiticity, gauge invariance, trace cyclicity) verified in unit tests.
4. **References:** Authoritative citations with DOI, arXiv ID, and citekeys matching `docs/references.bib` and the Obsidian Research Vault.

---

## 4. Integration with Code
* **In C++ Headers:** Public headers cite theory notes using `@see Theory: docs/theory/<note>.md#<section>`.
* **In Unit Tests:** Invariant tests cite equation numbers and verify analytical variations against brute-force recalculation.

