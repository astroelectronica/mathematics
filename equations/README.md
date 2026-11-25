# Equation Solver for GNU Octave

An interactive command-line program written for **GNU Octave** that solves equations and draws plots of the results.

The program asks for the **number of unknowns** and the **degree of the equation**, and then works in one of two modes:

| Unknowns | Degree | What the program does |
|---|---|---|
| any (1, 2, 3, ...) | 1 | Solves a **system of linear equations** `A*x = b` with several methods and many plot types |
| 1 | 2 or higher | Solves a **polynomial equation** `a_d*x^d + ... + a_1*x + a_0 = 0` and plots it |
| 2 or more | 2 or higher | Not supported (the program explains this and asks again) |

> A linear equation has degree 1, so the degree prompt decides between the linear solver and the polynomial solver.

## Requirements

- GNU Octave (a recent version is recommended).
- No extra packages are needed; the program only uses functions from Octave core.
- A graphical toolkit for the figures (the default Octave GUI/Qt or FLTK is fine).

## How to run

From a terminal:

```bash
octave linear_equation_solver.m
```

Or from the Octave prompt (recommended if the figure windows freeze while the program waits for input):

```octave
>> linear_equation_solver
```

## Linear mode (degree 1)

After choosing degree 1 the program asks for the **number of equations** and then, equation by equation, for the coefficient of each unknown and the constant term. The system can be square, overdetermined (more equations than unknowns) or underdetermined (fewer equations than unknowns).

### Menu

| Option | Description |
|---|---|
| 1 | Show the system in readable form, e.g. `2*x1 + 1*x2 = 5` |
| 2 | Analyze the system: ranks, determinant, condition number and type of solution |
| 3 | Solve the system (choose a method) |
| 4 | Verify the solution: shows `A*x`, `b` and the residual of every equation |
| 5 | Plot (choose the plot type) |
| 6 | Enter a new problem |
| 0 | Exit |

### Classification of the system

Option 2 uses the **Rouche-Capelli theorem**, comparing `rank(A)` and `rank([A|b])`:

| Condition | Result |
|---|---|
| `rank(A) < rank([A\|b])` | No solution (inconsistent) |
| `rank(A) = rank([A\|b]) = number of unknowns` | Exactly one solution |
| `rank(A) = rank([A\|b]) < number of unknowns` | Infinitely many solutions (`n - rank` free variables) |

### Solution methods

| Method | Applies to | Notes |
|---|---|---|
| 1. Backslash `A \ b` | Systems with a unique solution | Octave's built-in solver |
| 2. Gauss-Jordan (`rref`) | Any consistent system | Shows the reduced matrix. For infinite solutions it prints a particular solution and the direction vectors of the free variables (general solution `x = xp + t1*v1 + ...`) |
| 3. Cramer's rule | Square, `det(A) ~= 0` | Shows `det(A)` and every `det(A_i)` |
| 4. Inverse matrix | Square, `det(A) ~= 0` | `x = inv(A) * b` |
| 5. LU decomposition | Square, `det(A) ~= 0` | Shows `L`, `U`, `P` and solves with forward/backward substitution |
| 6. Pseudoinverse | Any system | `x = pinv(A) * b`: the least-squares solution if there is no exact solution, the minimum-norm solution if there are infinitely many |

If a method does not apply to the current system, the program says why and suggests an alternative.

### Plot types

| Option | Plot |
|---|---|
| 1 | **Geometric view**: curves for 1 unknown, lines for 2 unknowns, semi-transparent planes for 3 unknowns, with the computed solution as a black dot |
| 2 | Solution vector (bar chart, values written on the bars) |
| 3 | Solution vector (stem chart) |
| 4 | `A*x` compared with `b` (grouped bars) |
| 5 | Residual `A*x - b` per equation (bar chart) |
| 6 | Heat map of the augmented matrix `[A\|b]` |
| 7 | Smooth surface of `[A\|b]` (needs at least 2 equations) |
| 8 | 3D bars of `[A\|b]` (falls back to `stem3` if `bar3` is missing) |
| 9 | Singular values of `A` |
| 10 | Dashboard: six of the plots above in one figure |

If no solution has been computed yet, the plots use `pinv(A)*b` and the program prints a note. After each figure the program offers to save it as a PNG (for example `linear_system_geometry.png`).

## Polynomial mode (degree 2 or higher, one unknown)

The program asks for the coefficients from the highest power to the constant term (the leading coefficient cannot be zero).

| Option | Description |
|---|---|
| 1 | Show the equation, e.g. `1*x^2 - 3*x + 2 = 0` |
| 2 | Solve it with `roots`: real and complex roots, the discriminant for degree 2, and a check of `max |f(root)|` |
| 3 | Plot (choose the plot type) |
| 4 | Enter a new problem |
| 0 | Exit |

Plot types: the curve `f(x)` with its real roots, the roots in the complex plane, both side by side, and `f(x)` with its first two derivatives. The vertical range is chosen automatically so that the roots and the turning points stay visible.

## Example session

Linear system `2x + y = 5`, `x - y = 1` (solution `x1 = 2`, `x2 = 1`):

```
Number of unknowns: 2
Degree of the equation (1 = linear): 1
Number of equations: 2
  Equation 1: coefficients 2, 1 and constant 5
  Equation 2: coefficients 1, -1 and constant 1
Choose an option: 3        (solve)
Method: 3                  (Cramer's rule)
```

Polynomial `x^2 - 3x + 2 = 0` (roots `2` and `1`): enter `1` unknown, degree `2`, and the coefficients `1`, `-3`, `2`.

## Code structure

The file is an Octave **script** (it starts with `1;`) in which all functions are defined first and the main program is at the end.

| Section | Functions |
|---|---|
| Input helpers | `read_positive_int`, `read_number`, `offer_save` (every input is validated and asked again if invalid) |
| Linear systems | `read_system`, `print_system`, `classify_system`, `analyze_system`, `verify_solution` |
| Solution methods | `solve_system` |
| Linear plots | `draw_geometry`, `draw_solution_bar`, `draw_solution_stem`, `draw_ax_vs_b`, `draw_residual`, `draw_heatmap`, `draw_aug_surface`, `draw_aug_bars3d`, `draw_singular`, `plot_linear` |
| Linear mode | `linear_mode` (menu loop) |
| Polynomials | `poly_to_string`, `solve_polynomial`, `real_roots`, `poly_ymax`, `draw_poly_curve`, `draw_roots_plane`, `draw_poly_derivatives`, `plot_polynomial`, `polynomial_mode` |
| Main program | Asks for unknowns and degree and calls `linear_mode` or `polynomial_mode` until the user exits |

Design notes:

- Each `draw_*` function draws into the **current axes**, so the same function is used for a single figure and for a subplot of the dashboard.
- `linear_mode` and `polynomial_mode` return `true` when the user wants a new problem and `false` when the user exits.
- The input is read as text and converted with `str2double`, so invalid entries never crash the program.

## Limitations

- Equations of degree 2 or higher are only supported with one unknown. Nonlinear systems with several unknowns are not solved.
- Roots of polynomials are computed numerically; for repeated roots, tiny rounding errors can appear (a root such as `1e-6i` is treated as real when it is smaller than a tolerance of `1e-6`).
- The geometric view only exists for 1, 2 or 3 unknowns.
- Transparency of the planes (`FaceAlpha`) and legends of 3D surfaces depend on the graphics toolkit; if unsupported, the plot is still drawn without them.
- The program was written without being run on a real Octave installation, so if you find any compatibility problem with your Octave version, please report the error message.
