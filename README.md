# ARYA
ARYA is a CFD solver that is currently tested on a canonical lid-driven cavity flow (`sliding lid') problem on a uniform mesh using a finite-volume method with a first order upwind advection scheme. The sliding lid is shown in Figure 1. Momentum is diffused from the slip wall and eventually forms a large vortex that settles in the domain, and as the Reynolds number increases, secondary vortices form at the bottom corners. The solver uses a Semi-Implicit Method for Pressure Linked Equations (SIMPLE) approach. Equations are solved using a Tridiagonal Matrix Algorithm (TDMA).

<p align="center">
<img width="340" height="298" alt="sliding_lid drawio" src="https://github.com/user-attachments/assets/41033991-350b-496a-8e2c-dbfd4ee1b6c6" />
</p>

<p align="center">
  <b>Figure 1:</b> Sliding lid configuration.
</p>

## How to Use:
