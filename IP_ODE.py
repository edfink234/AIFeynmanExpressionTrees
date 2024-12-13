import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from itertools import product
import os

# Define the system of ODEs
def pendulum_derivatives(t, y, g, l, A, omega):
    theta, theta_dot = y
    dtheta_dt = theta_dot
    dtheta_dot_dt = - (g / l) * np.sin(theta) + (A * omega**2 / l) * np.cos(theta) * np.cos(omega * t)
    return [dtheta_dt, dtheta_dot_dt]

# Parameters
g = 9.81  # Acceleration due to gravity (m/s^2)
theta0 = np.pi / 6  # Initial angle (30 degrees)
theta_dot0 = 0.0    # Initial angular velocity
y0 = [theta0, theta_dot0]

# Time settings
t_start = 0.0
t_end = 20.0
num_points = 1000
t_eval = np.linspace(t_start, t_end, num_points)

# Parameter ranges
A_values = np.linspace(1, 10, 5)
omega_values = np.linspace(1, 10, 5)
l_values = np.linspace(1, 10, 5)

# Create output directory for figures
output_dir = "IP_figures"
os.makedirs(output_dir, exist_ok=True)

# Generate permutations of parameters
param_combinations = list(product(A_values, omega_values, l_values))

# Run simulations and save figures
time_image_paths = []  # Store paths of saved time images
phase_image_paths = []  # Store paths of saved phase images

for i, (A, omega, l) in enumerate(param_combinations):
    solution = solve_ivp(
        pendulum_derivatives,
        [t_start, t_end],
        y0,
        args=(g, l, A, omega),
        t_eval=t_eval
    )
    theta = solution.y[0]
    theta_dot = solution.y[1]
    
    print("iter =", i)
    # Plot theta vs time
    plt.figure(figsize=(12, 6))
    plt.plot(solution.t, theta, label=r'$\theta(t)$')
    plt.title(f'Pendulum Simulation ($A={A}, \omega={omega}, l={l}$)')
    plt.xlabel('Time (s)')
    plt.ylabel('Angle (rad)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    time_plot_path = os.path.join(output_dir, f'time_plot_{i}.png')
    plt.savefig(time_plot_path, dpi=5*96)
    plt.close()
    time_image_paths.append(time_plot_path)

    # Plot phase space (theta vs theta_dot) with vector field
    plt.figure(figsize=(8, 6))

    # Generate a grid for vector field
    
#    theta_grid, theta_dot_grid = np.meshgrid(theta[::20], theta_dot[::20])

    # Compute derivatives for the vector field
#    dtheta = theta_dot
#    dtheta_dot = - (g / l) * np.sin(theta) + (A * omega**2 / l) * np.cos(theta) * np.cos(omega * t_start)

    # Normalize vectors for better visualization
#    norm = np.sqrt(dtheta**2 + dtheta_dot**2)
#    dtheta = dtheta.copy() / norm
#    dtheta_dot = dtheta_dot.copy() / norm
    
    disperse_factor = 2
    # Plot vector field
#    plt.quiver(theta[::disperse_factor], theta_dot[::disperse_factor], dtheta[::disperse_factor], dtheta_dot[::disperse_factor], alpha = 0.2)
    
    dx = -A * omega * np.sin(omega * t_eval) + l * theta_dot * np.cos(theta)  # Compute dx/dt
    dy = -l * theta_dot * np.sin(theta)  # Compute dy/dt
    
    plt.quiver(theta[::disperse_factor], theta_dot[::disperse_factor], dx[::disperse_factor], dy[::disperse_factor], alpha = 0.2)


    # Plot trajectory
    plt.plot(theta, theta_dot, label='Trajectory')
    plt.title(f'Phase Space Plot with Vector Field ($A={A}, \omega={omega}, l={l}$)')
    plt.xlabel(r'$\theta$ (rad)')
    plt.ylabel(r'$\dot{\theta}$ (rad/s)')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    phase_plot_path = os.path.join(output_dir, f'phase_plot_{i}.png')
    plt.savefig(phase_plot_path, dpi=5*96)
    plt.close()
    phase_image_paths.append(phase_plot_path)
    
# Create HTML gallery
html_content = r"""
<html>
<head>
    <script type="text/javascript" defer
                src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML">
    </script>
    
    <style>
            body { font-family: Arial, sans-serif; }
            .gallery { display: flex; flex-wrap: wrap; }
            .gallery img { margin: 5px; width: calc(50% - 10px); } /* 5 images per row */
            h1 { text-align: center; }
            .equation { text-align: center; font-size: 1.2em; margin-bottom: 20px; }
    </style>
</head>
<body>
    <h1>Pendulum Simulation Gallery</h1>
    <div class="equation">
            $$\ddot{\theta} + \frac{g}{l} \sin (\theta) =  \frac{A \omega}{l} \left( \cos (\theta) \omega \cos(\omega t) \right), \; \theta_0 = \frac{\pi}{6}, \; \dot{\theta}_0 = 0$$
        </div>
    <div class="gallery">
"""

# Add images to the HTML
for path in time_image_paths:
    html_content += f'<img src="{os.path.basename(path)}" alt="Pendulum Simulation Figure">\n'
for path in phase_image_paths:
    html_content += f'<img src="{os.path.basename(path)}" alt="Pendulum Simulation Figure">\n'

html_content += """
    </div>
</body>
</html>
"""

# Save HTML file
html_path = os.path.join(output_dir, "gallery.html")
with open(html_path, "w") as f:
    f.write(html_content)

print(f"Gallery saved to {html_path}. Open this file in a browser to view the results.")

