import os
import glob
import json
import re

def extract_parts(filename):
  """
  Extracts parts of the form '{integer}_point_{integer}' from the given filename.

  Args:
    filename: The filename string.

  Returns:
    A list of extracted parts.
  """
  pattern = r"(\d+)_point_(\d+)"
  matches = re.findall(pattern, filename)
  parts = ["{}_point_{}".format(match[0], match[1]) for match in matches]
  return parts

def string_to_float(string_value):
    parts = string_value.split("_point_")
    if len(parts) != 2:
      return None  # Invalid format

    integer_part, decimal_part = parts
    return float(integer_part + "." + decimal_part)

def flt_to_str(flt):
    return str(flt).replace(".","_point_")

base_file = "trajectory_data_IC_2_point_2258_1_point_0_1_point_0_1_point_0_0_point_2_.mp4"
base_file_png = base_file.replace(".mp4",".png")
base_parts = [string_to_float(part) for part in extract_parts(base_file)]
#print(base_parts)

file_paths = glob.glob("*.*")
mp4_files = [f for f in file_paths if f.endswith(".mp4")]
png_files = [f for f in file_paths if f.endswith(".png")]
file_parts = []
file_float_parts = []

for file in mp4_files: #could be png_files as well, doesn't matter in this case.
    parts = extract_parts(file)
    file_parts.append(parts)
    file_float_parts.append([string_to_float(part) for part in parts])
#    print(file, file_parts[-1], file_float_parts[-1])
#for file, mp4_file in zip(file_float_parts, mp4_files):
#    print(file, mp4_file)

#print(f"len(file_float_parts) = {len(file_float_parts)}")
file_parts_display_1 = []
for parts in file_float_parts:
    if parts[-3:] == base_parts[-3:] and parts != base_parts:
        file_parts_display_1.append(parts)

file_parts_display_2 = []
for parts in file_float_parts:
    if [parts[1]] + parts[3:5] == [base_parts[1]] + base_parts[3:5] and parts != base_parts:
        file_parts_display_2.append(parts)

file_parts_display_3 = []
for parts in file_float_parts:
    if parts[1:3] + parts[4:] == base_parts[1:3] + base_parts[4:] and parts != base_parts:
        file_parts_display_3.append(parts)

file_parts_display_4 = []
for parts in file_float_parts:
    if parts[1:4] == base_parts[1:4] and parts != base_parts:
        file_parts_display_4.append(parts)

# Print the results (or use them as needed)
#print(f"file_parts_display_1 of length {len(file_parts_display_1)}:", *file_parts_display_1, sep = '\n')
#print(f"file_parts_display_2 of length {len(file_parts_display_2)}:", *file_parts_display_2, sep = '\n')
#print(f"file_parts_display_3 of length {len(file_parts_display_3)}:", *file_parts_display_3, sep = '\n')
#print(f"file_parts_display_4 of length {len(file_parts_display_4)}:", *file_parts_display_4, sep = '\n')
file_parts_display_1.sort(key = lambda x: x[1])
file_parts_display_2.sort(key = lambda x: x[2])
file_parts_display_3.sort(key = lambda x: x[3])
file_parts_display_4.sort(key = lambda x: x[4])
#print(f"file_parts_display_1 of length {len(file_parts_display_1)}:", *file_parts_display_1, sep = '\n')
#print(f"file_parts_display_2 of length {len(file_parts_display_2)}:", *file_parts_display_2, sep = '\n')
#print(f"file_parts_display_3 of length {len(file_parts_display_3)}:", *file_parts_display_3, sep = '\n')
#print(f"file_parts_display_4 of length {len(file_parts_display_4)}:", *file_parts_display_4, sep = '\n')
first_part = "trajectory_data_IC_"

file_parts_display_1_mp4 = [base_file]
file_parts_display_1_png = [base_file_png]
file_parts_display_2_mp4 = [base_file]
file_parts_display_2_png = [base_file_png]
file_parts_display_3_mp4 = [base_file]
file_parts_display_3_png = [base_file_png]
file_parts_display_4_mp4 = [base_file]
file_parts_display_4_png = [base_file_png]

for file_part in file_parts_display_1:
    temp_str = f"{first_part}{flt_to_str(file_part[0])}_{flt_to_str(file_part[1])}_{flt_to_str(file_part[2])}_{flt_to_str(file_part[3])}_{flt_to_str(file_part[4])}_"
    file_parts_display_1_mp4.append(f"{temp_str}.mp4")
    file_parts_display_1_png.append(f"{temp_str}.png")

for file_part in file_parts_display_2:
    temp_str = f"{first_part}{flt_to_str(file_part[0])}_{flt_to_str(file_part[1])}_{flt_to_str(file_part[2])}_{flt_to_str(file_part[3])}_{flt_to_str(file_part[4])}_"
    file_parts_display_2_mp4.append(f"{temp_str}.mp4")
    file_parts_display_2_png.append(f"{temp_str}.png")

for file_part in file_parts_display_3:
    temp_str = f"{first_part}{flt_to_str(file_part[0])}_{flt_to_str(file_part[1])}_{flt_to_str(file_part[2])}_{flt_to_str(file_part[3])}_{flt_to_str(file_part[4])}_"
    file_parts_display_3_mp4.append(f"{temp_str}.mp4")
    file_parts_display_3_png.append(f"{temp_str}.png")

for file_part in file_parts_display_4:
    temp_str = f"{first_part}{flt_to_str(file_part[0])}_{flt_to_str(file_part[1])}_{flt_to_str(file_part[2])}_{flt_to_str(file_part[3])}_{flt_to_str(file_part[4])}_"
    file_parts_display_4_mp4.append(f"{temp_str}.mp4")
    file_parts_display_4_png.append(f"{temp_str}.png")

titles = open("result_times.txt").readlines()
#print(titles)
display_1_titles = ["Base Case"] + [i.strip('\n') for i in titles[:10]]
display_2_titles = ["Base Case"] + [i.strip('\n') for i in titles[10:20]]
display_3_titles = ["Base Case"] + [i.strip('\n') for i in titles[20:30]]
display_4_titles = ["Base Case"] + [i.strip('\n') for i in titles[30:]]

print(*display_1_titles, sep='\n', end="\n\n")
print(*display_2_titles, sep='\n', end="\n\n")
print(*display_3_titles, sep='\n', end="\n\n")
print(*display_4_titles, sep='\n', end="\n\n")

display_1_title_tokens = [i.split() for i in display_1_titles[1:]]
display_2_title_tokens = [i.split() for i in display_2_titles[1:]]
display_3_title_tokens = [i.split() for i in display_3_titles[1:]]
display_4_title_tokens = [i.split() for i in display_4_titles[1:]]

display_1_titles = ["Base Case"]
display_2_titles = ["Base Case"]
display_3_titles = ["Base Case"]
display_4_titles = ["Base Case"]

print(*display_1_title_tokens, sep='\n', end="\n\n")
print(*display_2_title_tokens, sep='\n', end="\n\n")
print(*display_3_title_tokens, sep='\n', end="\n\n")
print(*display_4_title_tokens, sep='\n', end="\n\n")

for token_title in display_1_title_tokens:
    temp_str = f"Time spent learning {token_title[2]} = {token_title[6]} from {token_title[2]} = {token_title[4]} = {token_title[8]} seconds"
    if "not" in token_title:
        temp_str += r", \(\mathcal{L} = " f"{float(token_title[-1]):.4f}" r"\)"
    else:
        temp_str += r", \(\mathcal{L} < 1.3 \times 10^{-2}\)"
    display_1_titles.append(temp_str)

for token_title in display_2_title_tokens:
    temp_str = f"Time spent learning {token_title[2]} = {token_title[6]} from {token_title[2]} = {token_title[4]} = {token_title[8]} seconds"
    if "not" in token_title:
        temp_str += r", \(\mathcal{L} = " f"{float(token_title[-1]):.4f}" r"\)"
    else:
        temp_str += r", \(\mathcal{L} < 1.3 \times 10^{-2}\)"
        
    display_2_titles.append(temp_str)

for token_title in display_3_title_tokens:
    temp_str = f"Time spent learning {token_title[2]} = {token_title[6]} from {token_title[2]} = {token_title[4]} = {token_title[8]} seconds"
    if "not" in token_title:
        temp_str += r", \(\mathcal{L} = " f"{float(token_title[-1]):.4f}" r"\)"
    else:
        temp_str += r", \(\mathcal{L} < 1.3 \times 10^{-2}\)"
    display_3_titles.append(temp_str)

for token_title in display_4_title_tokens:
    temp_str = f"Time spent learning {token_title[2]} = {token_title[6]} from {token_title[2]} = {token_title[4]} = {token_title[8]} seconds"
    if "not" in token_title:
        temp_str += r", \(\mathcal{L} = " f"{float(token_title[-1]):.4f}" r"\)"
    else:
        temp_str += r", \(\mathcal{L} < 1.3 \times 10^{-2}\)"
    display_4_titles.append(temp_str)

#exit()
#print(*(file_parts_display_1_mp4+file_parts_display_1_png+file_parts_display_2_mp4+file_parts_display_2_png+file_parts_display_3_mp4+file_parts_display_3_png+file_parts_display_4_mp4+file_parts_display_4_png), sep='\n')
#[os.system(f"ls {i}") for i in (file_parts_display_1_mp4+file_parts_display_1_png+file_parts_display_2_mp4+file_parts_display_2_png+file_parts_display_3_mp4+file_parts_display_3_png+file_parts_display_4_mp4+file_parts_display_4_png)]

# Generate HTML
def generate_html(file_data):
    html_content = """<!DOCTYPE html>
<html lang="en">
<head>
    <script type="text/javascript" id="MathJax-script" async
        src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js">
    </script>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Generated Control Displays</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            display: flex;
            flex-direction: column;
            align-items: center;
            background-color: #f0f0f0;
        }
        .display {
            margin: 20px 0;
            display: flex;
            align-items: center;
            gap: 10px;
        }
        img, video {
            width: 45%;
            border: 1px solid #ccc;
            border-radius: 5px;
        }
        button {
            padding: 10px;
            background-color: #007bff;
            color: white;
            border: none;
            border-radius: 5px;
            cursor: pointer;
        }
        button:hover {
            background-color: #0056b3;
        }
        .title {
            font-size: 18px;
            font-weight: bold;
            text-align: center;
            margin-bottom: 10px;
        }
        .latex-equation {
            position: relative;
            display: inline-block;
        }

        .latex-equation.underline::after {
            content: "";
            position: absolute;
            left: 0;
            right: 0;
            bottom: 0;
            height: 1px; /* Thickness of the underline */
            background-color: currentColor; /* Matches the text color */
        }
        
        h3.latex-equation {
            margin-bottom: 0.2em; /* Add a small gap if needed */
        }

        div.latex-equation {
            margin-top: 0; /* Ensure there's no extra space above the div */
        }
    </style>
</head>
<body>
<h1>Control Problem Displays</h1>
"""r"""
<h2><u>Current Problem Formulation </u></h2>
<div class="latex-equation">
    $$\begin{align*}
    m \ddot{x} &= -\frac{d}{dx} \left(V_{\text{MT}}(x) + V_{\text{sech}}(x-\xi) \right) \\
    V_{\text{MT}}(x) &= \frac{1}{2} \Omega^2 x^2 \\
    V_{\text{sech}}(x-\xi) &= A \operatorname{sech}^{2}\left(b \left(x - {\xi}\right)\right)
    \end{align*}$$
</div>
<h3><u>Constraints</u></h3>
<div class="latex-equation">
    $$\begin{align*}
    x(0) &= x_0 \\
    \xi(0) &= 0 \\
    \dot{x}(0) &= 0 \\
    x(t \leq T) &= x^* = 0 \\
    \xi(t \leq T) &= x^* = 0 \\
    |\dot{x}(t \leq T)| &= 0\\
    T_{\text{th}} &= 10\text{ seconds} \\
    \end{align*}$$
</div>
<h3 class="latex-equation underline">Loss Function \(\mathcal{L}\)</h3>
<div class="latex-equation underline">
    $$\begin{align*}
    &\overrightarrow{\text{MSE}} \equiv \left(x^* - \vec{x}\right)\odot\left(x^* - \vec{x}\right) + \vec{v}\odot\vec{v} + \left(x^* - \vec{\xi}\right)\odot\left(x^* - \vec{\xi}\right) + \xi_{t=1}^2 \\
    &\text{MSE}_{\text{best}} = \min\left(\overrightarrow{\text{MSE}}\right) \\
    &t_{\text{best}} = \underset{t_i\, \in \,\vec{t}}{\text{argmin}} \left(\overrightarrow{\text{MSE}}\right) \\
    &\text{Smoothness Penalty} = \frac{1}{t_{\text{best}}}\sum_{t=1}^{t_{\text{best}}} \left(\xi_{t} - \xi_{t-1}\right)^2 \\
    &\text{Max}_{v} = \max\left(\left|\vec{v}_{[t=1:t_{\text{best}}]}\right|\right) \\
    &\text{Max}_{\xi} = \max\left(\left|\vec{\xi}_{[t=1:t_{\text{best}}]}\right|\right) \\[0.5cm]
    &\boxed{\mathcal{L} = \text{MSE}_{\text{best}} + \alpha \cdot \text{Smoothness Penalty} + \beta \cdot t_{\text{best}} + \gamma \cdot \text{Max}_{v} + \delta \cdot \text{Max}_{\xi}}
    \end{align*}$$
    <br> Where:
    $$\begin{align*}
    &\alpha = \beta = \gamma = \delta = 10^{-3} \\
    &\vec{x} = [x_1, x_2, ..., x_{T_{\text{Th}}}] \\
    &\vec{v} = [v_1, v_2, ..., v_{T_{\text{Th}}}] \\
    &\vec{\xi} = [\xi_1, \xi_2, ..., \xi_{T_{\text{Th}}}] 
    \end{align*}$$
</div>
<br>
<h3><u>Displays</u></h3>
"""
    for i, display in enumerate(file_data):
        first_png = display["png"][display["index"]]
        first_mp4 = display["mp4"][display["index"]]
        html_content += f"""
        <h4 id="title{i}" class="title latex-equation">{file_data[i]['titles'][file_data[i]['index']]}</h4> 
        <div class="display" id="display{i}">
            <button onclick="cycleDisplay({i}, -1)">&#8592;</button>
            <img id="img{i}" src="{first_png}" alt="Display {i + 1} Image">
            <video id="video{i}" src="{first_mp4}" controls></video>
            <button onclick="cycleDisplay({i}, 1)">&#8594;</button>
        </div>
        """

    html_content += r"""
<script>
    window.onload = function() {
        const displays = document.querySelectorAll('.display');
        displays.forEach(display => {
            const video = display.querySelector('video');
            const img = display.querySelector('img');
            
            // Ensure video metadata is loaded to access dimensions
            video.onloadedmetadata = () => {
                img.style.height = `${video.videoHeight / video.videoWidth * img.clientWidth}px`;
            };

            // If metadata is already loaded (e.g., cached), set height immediately
            if (video.readyState >= 1) {
                img.style.height = `${video.videoHeight / video.videoWidth * img.clientWidth}px`;
            }
        });
    };
    const fileData = """ + json.dumps(file_data) + r""";

    function cycleDisplay(displayIndex, direction) {
        const data = fileData[displayIndex];
        const currentIndex = data.index;
        const newIndex = (currentIndex + direction + data.png.length) % data.png.length;

        // Update image and video sources
        document.getElementById(`img${displayIndex}`).src = data.png[newIndex];
        document.getElementById(`video${displayIndex}`).src = data.mp4[newIndex];

        // Update the title
        document.getElementById(`title${displayIndex}`).innerText = data.titles[newIndex];
        // Trigger MathJax to re-render the LaTeX in the updated element
        MathJax.typesetPromise([document.getElementById(`title${displayIndex}`)]).catch((err) => console.error(err));
        // Update index
        fileData[displayIndex].index = newIndex;
    }
</script>
</body>
</html>
"""
    return html_content

# Prepare data for HTML generation
file_data = [
    {"png": file_parts_display_1_png, "mp4": file_parts_display_1_mp4, "index": 0, "titles": display_1_titles},
    {"png": file_parts_display_2_png, "mp4": file_parts_display_2_mp4, "index": 0, "titles": display_3_titles},
    {"png": file_parts_display_3_png, "mp4": file_parts_display_3_mp4, "index": 0, "titles": display_2_titles},
    {"png": file_parts_display_4_png, "mp4": file_parts_display_4_mp4, "index": 0, "titles": display_4_titles},
]

# Generate and save HTML
html_content = generate_html(file_data)
with open("control_displays.html", "w") as f:
    f.write(html_content)

print("HTML file generated as 'control_displays.html'")
