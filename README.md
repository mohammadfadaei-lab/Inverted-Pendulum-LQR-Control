# Inverted-Pendulum-LQR-Control
# Inverted Pendulum Control Using Linear Quadratic Regulator (LQR)

This repository provides a comprehensive MATLAB implementation of an **LQR (Linear Quadratic Regulator)** controller for the classical **cart–pole** (inverted pendulum) system. The code covers the entire design cycle: system modelling, controller synthesis, closed‑loop simulation, visualisation, sensitivity analysis, and a comparative study with a PID controller.

---

## Overview

The inverted pendulum is a canonical benchmark in control theory, characterised by inherent instability, underactuation, and nonlinear dynamics. This implementation utilises a linearised state‑space model around the upright equilibrium and applies an optimal LQR control law to stabilise the system. The software is structured to be modular, well‑commented, and easily adaptable for educational or research purposes.

---

## Key Features

- **Full state‑space modelling** of the cart–pendulum system, including friction and inertia terms.
- **LQR design** with adjustable weighting matrices \( \mathbf{Q} \) and \( R \) to trade off state regulation against control effort.
- **Time‑domain simulation** of the closed‑loop response under non‑zero initial conditions.
- **Comprehensive visualisation**:
  - State trajectories (cart position/velocity, pendulum angle/angular velocity).
  - Control signal (applied force).
  - Pole maps comparing open‑loop and closed‑loop eigenvalues.
- **Interactive animation** of the pendulum motion during simulation.
- **Sensitivity analysis** – investigates the effect of varying pendulum mass on closed‑loop performance.
- **Performance comparison** with a classical PID controller tuned for angle regulation.

---

## Prerequisites

- **MATLAB** (R2018a or later is recommended).
- **Control System Toolbox** – required for functions such as `lqr`, `ss`, `lsim`, `tf`, and `pid`.

---

## Repository Structure

- `inverted_pendulum_control.m` – the main MATLAB script containing all modelling, control design, simulation, and plotting routines.
- `README.md` – this file.
- `LICENSE` – the MIT license under which this project is released.

---

## Usage Instructions

1. Clone or download the repository to your local machine.
2. Open MATLAB and navigate to the project folder.
3. Open the file `inverted_pendulum_control.m` in the MATLAB Editor.
4. Run the script (click the green **Run** button or type the file name in the Command Window).
5. The script will execute all sections sequentially and produce:
   - Printed diagnostic information in the Command Window.
   - Several figure windows showing simulation results, sensitivity analysis, and controller comparisons.
   - An animated visualisation of the stabilised pendulum.

---

## LQR Tuning and Customisation

The LQR cost function is defined as

\[
J = \int_0^\infty \left( \mathbf{x}^\top \mathbf{Q} \mathbf{x} + u^\top R u \right) dt,
\]

with the control law \( u = -\mathbf{K} \mathbf{x} \). The weighting matrices are set in Section 4 of the script:

```matlab
Q = [100,   0,   0,   0;    % cart position
       0,   1,   0,   0;    % cart velocity
       0,   0, 500,   0;    % pendulum angle
       0,   0,   0,  50];   % angular velocity

R = 1;                      % control effort penalty
```

- Increasing the diagonal entries of \( \mathbf{Q} \) (especially those corresponding to the cart position and pendulum angle) yields faster transient response but demands higher control forces.
- Increasing \( R \) reduces the control effort, resulting in more conservative (slower) regulation.

Users are encouraged to adjust these matrices to explore the trade‑off between performance and energy consumption.

---

## Simulation Results and Outputs

Upon execution, the script generates the following graphical outputs:

1. **State trajectories** (cart displacement, cart velocity, pendulum angle, angular velocity) for both the uncontrolled and LQR‑controlled systems.
2. **Control force** applied to the cart over the simulation horizon.
3. **Pole placement comparison** – open‑loop poles (unstable) versus closed‑loop poles (stable, with desirable damping).
4. **Sensitivity analysis** – response of the LQR‑controlled system for several values of pendulum mass, showing robustness to parameter variations.
5. **LQR vs PID comparison** – superimposed plots of angle regulation and control effort, illustrating the advantages of optimal control.
6. **Real‑time animation** – a moving cart and pendulum visualisation that updates during the simulation (frame rate is optimised for performance).

All simulation data are saved to `inverted_pendulum_results.mat` for post‑processing.

---

## License

This project is distributed under the **MIT License** – you are free to use, modify, and distribute it, provided that the original copyright notice is retained. See the `LICENSE` file for full details.

---

## Citation and Contributions

If you find this implementation useful for your research or teaching, please consider citing it by referencing the repository URL:

> *Inverted Pendulum Control Using LQR*, GitHub repository, [https://github.com/[YourUsername]/[YourRepoName]](https://github.com/mohammadfadaei-lab/Inverted-Pendulum-LQR-Control).

Contributions in the form of bug reports, suggestions, or pull requests are warmly welcomed. Please open an issue or submit a pull request via the GitHub interface.

---

**Developed with ❤️ and MATLAB.**
