Parametric Minimization Integration Test

Tests simultaneous minimization of Crick parametric DOFs (r0, omega0,
delta_omega0, etc.) alongside standard atom-tree DOFs (chi angles).

Creates a 3-helix bundle using MakeBundle with a realistic ELLKAIA heptad
repeat sequence (21 residues per helix, 63 total), then minimizes using
MinMover with the new parametric DOF support.

Two runs are performed:
  1. Parametric DOFs only (chi angles frozen)
  2. Parametric DOFs + chi angles simultaneously

Both runs should converge to lower energy than the starting structure.
The parametric+chi run should achieve at least as low energy as the
parametric-only run.
