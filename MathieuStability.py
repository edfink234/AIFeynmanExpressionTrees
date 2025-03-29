import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Parameters
T = 200        # maximum integration time
delta = 1      # maximum allowed deviation for stability
phi0 = 0.0     # initial phi value
dphi0_dt = 1e-4  # initial dphi/dt

# Parameter grid for k and m_tilde
k_values = np.arange(-1, 1.8, 0.01)
m_tilde_values = np.arange(0, 2, 0.01)
K, M = np.meshgrid(k_values, m_tilde_values, indexing='ij')
N = K.size
k_flat = K.flatten()
m_flat = M.flatten()

# Build the initial condition for all trajectories: shape (2*N,)
# Each trajectory has initial [phi0, dphi0_dt]
y0 = np.zeros(2 * N)
y0[0::2] = phi0      # phi initial conditions
y0[1::2] = dphi0_dt   # dphi/dt initial conditions

def vectorized_mathieu(t, y):
    """
    y is a vector of shape (2*N, ) where:
      y[0::2] are the phi_i values,
      y[1::2] are the dphi_i/dt values,
    and k_flat and m_flat are the parameter arrays for each trajectory.
    """
    # Extract phi and dphi/dt for all trajectories:
    phi = y[0::2]
    dphi_dt = y[1::2]
    
    # Compute second derivative for each trajectory:
    d2phi_dt2 = (-k_flat + m_flat * np.cos(t)) * phi
    
    # Build the derivative vector:
    dydt = np.empty_like(y)
    dydt[0::2] = dphi_dt
    dydt[1::2] = d2phi_dt2
    return dydt

# Integrate the full vectorized system
sol = solve_ivp(vectorized_mathieu, [0, T], y0, method='DOP853', dense_output=True)

# Evaluate the solution at a grid of times
t_eval = np.linspace(0, T, 1000)
Y = sol.sol(t_eval)  # shape: (2*N, len(t_eval))
phi_all = Y[0::2, :]  # shape: (N, len(t_eval)) for all phi trajectories

# For each trajectory, compute the maximum deviation from phi0
max_deviation = np.max(np.abs(phi_all - phi0), axis=1)

# Identify stable/unstable trajectories
stable_mask = max_deviation < delta
stable_points = np.column_stack((k_flat[stable_mask], m_flat[stable_mask]))
unstable_points = np.column_stack((k_flat[~stable_mask], m_flat[~stable_mask]))

# Plot the stability regions
plt.figure(figsize=(10, 6))
if stable_points.size > 0:
    plt.scatter(stable_points[:, 0], stable_points[:, 1], color='green', s=1, label='Stable')
if unstable_points.size > 0:
    plt.scatter(unstable_points[:, 0], unstable_points[:, 1], color='red', s=1, label='Unstable')
plt.xlabel(r'$k$')
plt.ylabel(r'$\tilde{m}$')
plt.title(rf'$\mathrm{{Max}}\left(|\phi-\phi_0|\right) < {delta}, T={T}, \dot{{\phi}}_0 = {dphi0_dt:.2e}$')
plt.legend()
plt.grid(True)
plt.show()
