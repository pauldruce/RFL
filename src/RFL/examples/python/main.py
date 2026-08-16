import rfl
import numpy as np

def main():
    print("========================================")
    print(" RFL Python Binding Example")
    print("========================================\n")
    
    # Configure safety limits
    default_max = rfl.get_max_clifford_mode()
    print(f"Default Max Clifford Mode (p+q): {default_max}")
    
    new_max = 20
    print(f"Overriding Max Clifford Mode to: {new_max} (allows larger matrices)\n")
    rfl.set_max_clifford_mode(new_max)

    # Initialize a Dirac Operator
    p, q = 1, 3
    dim = 10
    print(f"Initializing Dirac Operator with p={p}, q={q}, dim={dim}...")
    dirac = rfl.DiracOperator(p, q, dim)

    print(f"-> Dirac Operator Type: {dirac.get_type()}")
    print(f"-> Matrix Dimension: {dirac.get_matrix_dimension()}\n")

    # Setup the Metropolis Algorithm
    g_2 = -1.0
    g_4 = 1.0
    scale = 1.0
    num_steps = 100
    seed = 42

    print(f"Running Metropolis Algorithm (g_2={g_2}, g_4={g_4}, steps={num_steps})...")
    metropolis = rfl.Metropolis(g_2, g_4, scale, num_steps, seed)
    acceptance_rate = metropolis.update_dirac(dirac)
    print(f"-> Acceptance Rate: {acceptance_rate * 100:.2f}%\n")

    # Analyze Eigenvalues
    print("Computing eigenvalues...")
    eigenvals = dirac.get_eigenvalues()
    print(f"-> Found {len(eigenvals)} eigenvalues.")
    print(f"-> Min eigenvalue: {np.min(eigenvals):.4f}")
    print(f"-> Max eigenvalue: {np.max(eigenvals):.4f}")
    print(f"-> Mean eigenvalue: {np.mean(eigenvals):.4f}")
    print("\nSimulation complete!")

if __name__ == "__main__":
    main()
