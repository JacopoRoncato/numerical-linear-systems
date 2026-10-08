[README (1).md](https://github.com/user-attachments/files/33211024/README.1.md)

# Numerical Linear Systems — MATLAB

**Numerical methods for solving a resistor-network linear system**, implemented in MATLAB as part of the *Fondamenti di Calcolo Numerico* (Foundations of Numerical Computing) course at Politecnico di Milano (2026).

This project formulates a nine-unknown linear system from a resistor network and compares direct and iterative solvers, focusing on convergence, numerical error, and conditioning.

## Methods

**Direct methods**
- LU factorization using MATLAB's `lu`
- Forward substitution (`fwsub.m`) and backward substitution (`bksub.m`)
- Comparison against MATLAB's backslash solver (`A\b`)

**Iterative methods**
- Jacobi iteration (`jacobi.m`)
- Gradient method / steepest descent (`graddyn.m`)
- Diagonally preconditioned gradient method (`gradprec.m`)
- Conjugate gradient (`gc.m`)

The script evaluates relative solution errors against the backslash reference, iteration counts, normalized residual histories, the Jacobi iteration matrix's spectral radius, and matrix condition numbers.

## Project structure

```text
numerical-linear-systems/
├── README.md
├── .gitignore
├── Progetto_1.m
├── fwsub.m
├── bksub.m
├── jacobi.m
├── graddyn.m
├── gradprec.m
└── gc.m
```

## Requirements

- MATLAB with standard numerical linear algebra and plotting functions
- No third-party packages are used in the supplied source files

## Run

1. Download or clone the repository.
2. Open MATLAB and set the **Current Folder** to the repository directory.
3. Run the main script:

```matlab
Progetto_1
```

Keep all six function files in the same folder as `Progetto_1.m`, or ensure they are on the MATLAB path.

## Numerical experiment

The script constructs a **9 × 9 conductance matrix** from a resistor network with `Vin = 15` and `Rin = 60`, then solves `A*x = b`.

The iterative solvers start from a zero vector, with a normalized-residual tolerance of `1e-11` and an iteration limit of `1000`.

The script prints numerical comparisons and opens a semilogarithmic plot of normalized residual versus iteration for the gradient, preconditioned gradient, and conjugate-gradient implementations.

### Convergence plot

After running the main script, the plot can be saved from the MATLAB figure window. To show it directly here on GitHub, save it as `results/convergence.png` and uncomment the line below:

<!-- ![Normalized residual convergence](results/convergence.png) -->

## Notes and limitations

- This is a course assignment and a **small educational numerical experiment**, not a general-purpose linear-system solver library.
- The six solver functions are retained from the supplied project implementations; the main script has been corrected for undefined variables, spectral-radius calculation, and convergence plotting.
- Execution and numerical equivalence have **not yet been verified in a MATLAB runtime** in the preparation of this repository. Results may depend on numerical precision and the implementations provided.
- The published script is intended for study and reproducibility, not for production electrical-network analysis.

## Academic context

**Course:** Fondamenti di Calcolo Numerico (2026)  
**Institution:** Politecnico di Milano  
**Language:** MATLAB
