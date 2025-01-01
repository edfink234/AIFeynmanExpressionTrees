import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import os
import pandas as pd
import glob
import json

# Constants
g = 9.8  # Gravity acceleration, m/s^2
t_span = (0, 2)  # Time range
t_points = 1000  # Number of time points
t_eval = np.linspace(t_span[0], t_span[1], t_points)  # Evaluation times

# Create output directory for figures
output_dir = "IP_COMPARISON_FIGURES"
os.makedirs(output_dir, exist_ok=True)

# Parameter ranges
A_values = np.linspace(1, 10, 6)
omega_values = np.linspace(1, 10, 6)
l_values = np.linspace(1, 10, 6)
theta_0_values = np.linspace(np.pi - 0.125, np.pi + 0.125, 6)
# Plot results with LaTeX-style labels
#labels = (r"$\ddot{\theta} + \frac{g}{l} \sin (\theta) =  \frac{ A \omega^2}{l} \left( \cos (\theta) \cos(\omega t) \right)$", r"$\ddot{\theta} + \frac{g}{l}(\pi - \theta) =  -\frac{ A \omega^2}{l}\cos(\omega t)$", r"$\ddot{\theta} + \frac{g}{l}\left((\pi - \theta) + \frac{(\theta - \pi)^3}{6}\right) =  \frac{ A \omega^2}{l}\left(-1 + \frac{(x-\pi)^2}{2}\right)\cos(\omega t)$")
labels = ('Original', 'First-Order', 'Second Order')

# Function to convert second-order ODE to a first-order system
def equation_system(t, y, params, equation_type):
    theta, theta_dot = y
    A, omega, l = params
    
    if equation_type == 1:  # Original equation
        theta_ddot = -(g / l) * np.sin(theta) + (A * omega**2 / l) * np.cos(theta) * np.cos(omega * t)
    elif equation_type == 2:  # First-order Taylor expansion
        theta_ddot = -(g / l) * (np.pi - theta) - (A * omega**2 / l) * np.cos(omega * t)
    elif equation_type == 3:  # Second-order Taylor expansion
        theta_ddot = -(g / l) * ((np.pi - theta) + (theta - np.pi)**3 / 6) + \
                     (A * omega**2 / l) * ((-1 + (theta - np.pi)**2 / 2) * np.cos(omega * t))
    return [theta_dot, theta_ddot]

# Create parameter mesh grid
parameter_grid = np.array(np.meshgrid(A_values, omega_values, l_values, theta_0_values)).T.reshape(-1, 4)

# Simulate for all parameter sets
for idx, (A, omega, l, theta_0) in enumerate(parameter_grid):
    output_svg = os.path.join(output_dir, f'A_{A}_omega{omega}_l_{l}_theta_0_{theta_0}.svg')
    output_csv = os.path.join(output_dir, f'A_{A}_omega{omega}_l_{l}_theta_0_{theta_0}.csv')
    
    # Skip if SVG already exists
    if os.path.exists(output_svg):
        print(f"Skipping simulation {idx+1}/{len(parameter_grid)}: Results already exist.")
        continue

    print(f"Simulating set {idx+1}/{len(parameter_grid)}: A={A}, omega={omega}, l={l}, theta_0={theta_0}")
    
    # Initial conditions
    y0 = [theta_0, 0]  # Initial [theta, theta_dot]

    # Solve equations
    solutions = []
    for eq_type in range(1, 4):
        sol = solve_ivp(equation_system, t_span, y0, t_eval=t_eval, args=([A, omega, l], eq_type))
        solutions.append(sol.y[0])  # Store theta(t)

    # Save solutions to CSV
    data = np.column_stack((t_eval, np.array(solutions).T))
    df = pd.DataFrame(data, columns=['t_eval', 'Equation_1', 'Equation_2', 'Equation_3'])
    df.to_csv(output_csv, index=False)

    # Plot results
    plt.figure(figsize=(10, 6))
    for theta_t, label in zip(solutions, labels):
        plt.plot(t_eval, theta_t, label=label)
    
    plt.title(rf'Simulation: A={A}, $\omega$={omega}, l={l}, $\theta_0$={theta_0}')
    plt.xlabel('Time (t)')
    plt.ylabel(r'$\theta(t)$')
    plt.legend()
    plt.grid()
    plt.savefig(output_svg, format='svg')  # Save the plot as SVG
    plt.close()

# Get all available SVG files in the directory
svg_files = glob.glob(os.path.join(output_dir, "*.svg"))
print(len(svg_files))

# Parse filenames to extract parameter values
images_data = []
for svg_file in svg_files:
    filename = os.path.basename(svg_file)
    try:
        # Extract parameters from the filename (e.g., A_0.1_omega0.1_l_0.1_theta_0_2.14.svg)
        parts = filename.replace(".svg", "").split("_")
        A = float(parts[1])
        omega = float(parts[2].replace('omega',''))
        l = float(parts[4])
        theta_0 = float(parts[7])
        images_data.append({"A": A, "omega": omega, "l": l, "theta_0": theta_0, "file":  os.path.join(output_dir, filename)})
        print(images_data[-1])
    except (IndexError, ValueError):
        print(f"Skipping invalid filename: {filename}")

# Prepare JavaScript variable from the images data
images_json = json.dumps(images_data)

# HTML template with left/right arrows and default values as minimums
html_template = r"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Simulation Viewer</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        .control-container { margin-bottom: 20px; }
        .control-label { margin-right: 10px; }
        .arrow { cursor: pointer; font-size: 18px; padding: 5px; }
        .slider { margin-left: 10px; width: 200px; }
        img { max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 4px; padding: 5px; }
    </style>
</head>
<body>
    <h1>Simulation Viewer</h1>

    <div class="control-container">
        <label class="control-label">A:</label>
        <span class="arrow" onclick="changeValue('A', -1)">←</span>
        <span id="A-value">0.0</span>
        <span class="arrow" onclick="changeValue('A', 1)">→</span>
        <input id="A-slider" type="range" min="0" max="0" step="1" class="slider" oninput="updateSliderValue('A', this.value)">
    </div>
    <div class="control-container">
        <label class="control-label">ω:</label>
        <span class="arrow" onclick="changeValue('omega', -1)">←</span>
        <span id="omega-value">0.0</span>
        <span class="arrow" onclick="changeValue('omega', 1)">→</span>
        <input id="omega-slider" type="range" min="0" max="0" step="1" class="slider" oninput="updateSliderValue('omega', this.value)">
    </div>
    <div class="control-container">
        <label class="control-label">l:</label>
        <span class="arrow" onclick="changeValue('l', -1)">←</span>
        <span id="l-value">0.0</span>
        <span class="arrow" onclick="changeValue('l', 1)">→</span>
        <input id="l-slider" type="range" min="0" max="0" step="1" class="slider" oninput="updateSliderValue('l', this.value)">
    </div>
    <div class="control-container">
        <label class="control-label">θ₀:</label>
        <span class="arrow" onclick="changeValue('theta_0', -1)">←</span>
        <span id="theta_0-value">0.0</span>
        <span class="arrow" onclick="changeValue('theta_0', 1)">→</span>
        <input id="theta_0-slider" type="range" min="0" max="0" step="1" class="slider" oninput="updateSliderValue('theta_0', this.value)">
    </div>

    <div>
        <img id="simulation-image" src="" alt="Simulation Image">
    </div>

    <script>
        // Available images data
        const images = {images_json};

        // Get unique values for each parameter
        const uniqueValues = {
            A: [...new Set(images.map(img => img.A))].sort((a, b) => a - b),
            omega: [...new Set(images.map(img => img.omega))].sort((a, b) => a - b),
            l: [...new Set(images.map(img => img.l))].sort((a, b) => a - b),
            theta_0: [...new Set(images.map(img => img.theta_0))].sort((a, b) => a - b),
        };

        // Initialize current values to the minimums
        const currentValues = {
            A: uniqueValues.A[0],
            omega: uniqueValues.omega[0],
            l: uniqueValues.l[0],
            theta_0: uniqueValues.theta_0[0],
        };

        // Initialize sliders
        function initSliders() {
            for (const param in uniqueValues) {
                const slider = document.getElementById(`${param}-slider`);
                const values = uniqueValues[param];
                slider.min = 0;
                slider.max = values.length - 1;
                slider.value = values.indexOf(currentValues[param]);
            }
        }

        // Update displayed image
        function updateImage() {
            const matchingImage = images.find(img =>
                Math.abs(img.A - currentValues.A) < 0.01 &&
                Math.abs(img.omega - currentValues.omega) < 0.01 &&
                Math.abs(img.l - currentValues.l) < 0.01 &&
                Math.abs(img.theta_0 - currentValues.theta_0) < 0.01
            );

            const imageElement = document.getElementById("simulation-image");
            if (matchingImage) {
                imageElement.src = matchingImage.file;
                imageElement.alt = `Simulation Image: A=${matchingImage.A}, omega=${matchingImage.omega}, l=${matchingImage.l}, theta_0=${matchingImage.theta_0}`;
            } else {
                imageElement.src = "";
                imageElement.alt = "No matching image found";
            }

            // Update displayed parameter values
            document.getElementById("A-value").textContent = currentValues.A.toFixed(2);
            document.getElementById("omega-value").textContent = currentValues.omega.toFixed(2);
            document.getElementById("l-value").textContent = currentValues.l.toFixed(2);
            document.getElementById("theta_0-value").textContent = currentValues.theta_0.toFixed(2);
        }

        // Change value for a parameter
        function changeValue(param, direction) {
            const values = uniqueValues[param];
            const currentIndex = values.indexOf(currentValues[param]);
            const newIndex = currentIndex + direction;

            if (newIndex >= 0 && newIndex < values.length) {
                currentValues[param] = values[newIndex];
                document.getElementById(`${param}-slider`).value = newIndex;
                updateImage();
            }
        }

        // Update parameter value from slider
        function updateSliderValue(param, sliderIndex) {
            const values = uniqueValues[param];
            currentValues[param] = values[sliderIndex];
            updateImage();
        }

        // Initialize display
        initSliders();
        updateImage();
    </script>
</body>
</html>
""".replace('{images_json}', f"{images_json}")

# Save HTML to file
output_html = "ip_comparison.html"
with open(output_html, "w") as f:
    f.write(html_template)

print(f"HTML page generated: {output_html}")
