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
base_file_1_point_4 = "trajectory_data_IC_2_point_3687_1_point_4_1_point_0_1_point_0_0_point_2_.mp4"
base_file_1_point_4_png = base_file_1_point_4.replace(".mp4",".png")
base_file_1_point_1 = "trajectory_data_IC_2_point_266292_1_point_1_1_point_0_1_point_0_0_point_2_.mp4"
base_file_1_point_1_png = base_file_1_point_1.replace(".mp4",".png")
base_parts = [string_to_float(part) for part in extract_parts(base_file)]
#print(base_parts)

file_paths = glob.glob("*.*")
is_normal = lambda f: 'random_initialized' not in f and 'Variational' not in f
mp4_files = [f for f in file_paths if f.endswith(".mp4") and is_normal(f)]
png_files = [f for f in file_paths if f.endswith(".png") and is_normal(f)]
print(len(mp4_files))
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
file_parts_display_2 = []
file_parts_display_3 = []
file_parts_display_4 = []
file_parts_display_5 = []
file_parts_display_6 = []
file_parts_display_7 = []
file_parts_display_8 = ['trajectory_data_IC_2_point_266292_1_point_1_1_point_0_1_point_0_0_point_2_.png', 'trajectory_data_IC_2_point_273769_1_point_099_0_point_99_1_point_0_0_point_2015_.png', 'trajectory_data_IC_2_point_281359_1_point_098_0_point_98_1_point_0_0_point_203_.png', 'trajectory_data_IC_2_point_289064_1_point_097_0_point_97_1_point_0_0_point_2045_.png', 'trajectory_data_IC_2_point_296884_1_point_096_0_point_96_1_point_0_0_point_206_.png', 'trajectory_data_IC_2_point_30482_1_point_095_0_point_95_1_point_0_0_point_2075_.png', 'trajectory_data_IC_2_point_312874_1_point_094_0_point_94_1_point_0_0_point_209_.png', 'trajectory_data_IC_2_point_321044_1_point_093_0_point_93_1_point_0_0_point_2105_.png', 'trajectory_data_IC_2_point_329333_1_point_092_0_point_92_1_point_0_0_point_212_.png', 'trajectory_data_IC_2_point_337739_1_point_091_0_point_91_1_point_0_0_point_2135_.png', 'trajectory_data_IC_2_point_346265_1_point_09_0_point_9_1_point_0_0_point_215_.png', 'trajectory_data_IC_2_point_354909_1_point_089_0_point_89_1_point_0_0_point_2165_.png', 'trajectory_data_IC_2_point_363673_1_point_088_0_point_88_1_point_0_0_point_218_.png', 'trajectory_data_IC_2_point_372557_1_point_087_0_point_87_1_point_0_0_point_2195_.png', 'trajectory_data_IC_2_point_38156_1_point_086_0_point_86_1_point_0_0_point_221_.png', 'trajectory_data_IC_2_point_390682_1_point_085_0_point_85_1_point_0_0_point_2225_.png', 'trajectory_data_IC_2_point_399923_1_point_084_0_point_84_1_point_0_0_point_224_.png', 'trajectory_data_IC_2_point_409283_1_point_083_0_point_83_1_point_0_0_point_2255_.png', 'trajectory_data_IC_2_point_418761_1_point_082_0_point_82_1_point_0_0_point_227_.png', 'trajectory_data_IC_2_point_428356_1_point_081_0_point_81_1_point_0_0_point_2285_.png', 'trajectory_data_IC_2_point_438068_1_point_08_0_point_8_1_point_0_0_point_23_.png']
file_parts_display_8 = [extract_parts(filename) for filename in file_parts_display_8]
file_parts_display_8 = [[string_to_float(part) for part in file] for file in file_parts_display_8]
file_parts_display_9 = ['trajectory_data_IC_2_point_438068_1_point_08_0_point_8_1_point_0_0_point_23_.png', 'trajectory_data_IC_2_point_434087_1_point_0_1_point_0_1_point_0_0_point_2_Variational_.png']
file_parts_display_9_mp4 = [i.replace('.png', '.mp4') for i in file_parts_display_9]
file_parts_display_9_png = file_parts_display_9.copy()
file_parts_display_9 = [extract_parts(filename) for filename in file_parts_display_9]
file_parts_display_9 = [[string_to_float(part) for part in file] for file in file_parts_display_9]
#print(*file_parts_display_8, sep='\n')

for parts in file_float_parts:
    if parts[-3:] == base_parts[-3:] and parts != base_parts and not (1.4 < parts[1] < 1.5):
        file_parts_display_1.append(parts)
    elif [parts[1]] + parts[3:5] == [base_parts[1]] + base_parts[3:5] and parts != base_parts:
        file_parts_display_2.append(parts)
    elif parts[1:3] + parts[4:] == base_parts[1:3] + base_parts[4:] and parts != base_parts:
        file_parts_display_3.append(parts)
    elif parts[1:4] == base_parts[1:4] and parts != base_parts:
        file_parts_display_4.append(parts)
    elif parts[-3:] == base_parts[-3:] and parts != base_parts and (1.4 < parts[1] < 1.5):
        file_parts_display_5.append(parts)
    elif parts[-2:] == base_parts[-2:] and parts != base_parts and (0.66 <= parts[1] <= 0.99) and (0.75 <= parts[2] <= 0.99):
        file_parts_display_6.append(parts)
    elif parts[-2:] == base_parts[-2:] and parts != base_parts and (0.66 <= parts[1] <= 0.68) and (0.67 <= parts[2] <= 0.74):
        file_parts_display_7.append(parts)

#all_variables = globals()
#print(sum([len(all_variables[i]) for i in all_variables if "file_parts_display_" in i]))
#exit()
# Print the results (or use them as needed)
file_parts_display_1 = [i for i in file_parts_display_1 if i != [2.266292, 1.1, 1.0, 1.0, 0.2]]
#print(f"file_parts_display_2 of length {len(file_parts_display_2)}:", *file_parts_display_2, sep = '\n')
#print(f"file_parts_display_3 of length {len(file_parts_display_3)}:", *file_parts_display_3, sep = '\n')
#print(f"file_parts_display_4 of length {len(file_parts_display_4)}:", *file_parts_display_4, sep = '\n')
#print(f"file_parts_display_5 of length {len(file_parts_display_5)}:", *file_parts_display_5, sep = '\n')
#print(f"file_parts_display_6 of length {len(file_parts_display_6)}:", *file_parts_display_6, sep = '\n')
#print(f"file_parts_display_7 of length {len(file_parts_display_7)}:", *file_parts_display_7, sep = '\n')

file_parts_display_1.sort(key = lambda x: x[1])
file_parts_display_2.sort(key = lambda x: x[2])
file_parts_display_3.sort(key = lambda x: x[3])
file_parts_display_4.sort(key = lambda x: x[4])
file_parts_display_5.sort(key = lambda x: x[1])
file_parts_display_6.sort(key = lambda x: x[1] + x[2])
#print("\n"*10)
#print(*list(zip(range(1, len(file_parts_display_6)+1),file_parts_display_6)), sep='\n')
#exit()
file_parts_display_6[-4], file_parts_display_6[-5] = file_parts_display_6[-5], file_parts_display_6[-4]
file_parts_display_6[-17], file_parts_display_6[-18] = file_parts_display_6[-18], file_parts_display_6[-17]
file_parts_display_6[-30], file_parts_display_6[-31] = file_parts_display_6[-31], file_parts_display_6[-30]
file_parts_display_6[-36], file_parts_display_6[-37] = file_parts_display_6[-37], file_parts_display_6[-36]
file_parts_display_7.sort(key = lambda x: x[1] + x[2])
#swapping last 2 files because for file_parts_display_7 it's not perfect
file_parts_display_7[0], file_parts_display_7[1] = file_parts_display_7[1], file_parts_display_7[0]
#print(*file_parts_display_7, sep='\n')
#exit()
#print(f"file_parts_display_7 of length {len(file_parts_display_7)}:", *file_parts_display_7, sep = '\n')
#exit()

first_part = "trajectory_data_IC_"

file_parts_display_1_mp4 = [base_file]
file_parts_display_1_png = [base_file_png]
file_parts_display_2_mp4 = [base_file]
file_parts_display_2_png = [base_file_png]
file_parts_display_3_mp4 = [base_file]
file_parts_display_3_png = [base_file_png]
file_parts_display_4_mp4 = [base_file]
file_parts_display_4_png = [base_file_png]
file_parts_display_5_mp4 = [base_file_1_point_4]
file_parts_display_5_png = [base_file_1_point_4_png]
file_parts_display_6_mp4 = [base_file]
file_parts_display_6_png = [base_file_png]
file_parts_display_7_mp4 = [base_file]
file_parts_display_7_png = [base_file_png]
file_parts_display_8_mp4 = []
file_parts_display_8_png = []

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

for file_part in file_parts_display_5:
    temp_str = f"{first_part}{flt_to_str(file_part[0])}_{flt_to_str(file_part[1])}_{flt_to_str(file_part[2])}_{flt_to_str(file_part[3])}_{flt_to_str(file_part[4])}_"
    file_parts_display_5_mp4.append(f"{temp_str}.mp4")
    file_parts_display_5_png.append(f"{temp_str}.png")

for file_part in file_parts_display_6[::-1]:
    temp_str = f"{first_part}{flt_to_str(file_part[0])}_{flt_to_str(file_part[1])}_{flt_to_str(file_part[2])}_{flt_to_str(file_part[3])}_{flt_to_str(file_part[4])}_"
    file_parts_display_6_mp4.append(f"{temp_str}.mp4")
    file_parts_display_6_png.append(f"{temp_str}.png")

for file_part in file_parts_display_7[::-1]:
    temp_str = f"{first_part}{flt_to_str(file_part[0])}_{flt_to_str(file_part[1])}_{flt_to_str(file_part[2])}_{flt_to_str(file_part[3])}_{flt_to_str(file_part[4])}_"
    file_parts_display_7_mp4.append(f"{temp_str}.mp4")
    file_parts_display_7_png.append(f"{temp_str}.png")

for file_part in file_parts_display_8:
    temp_str = f"{first_part}{flt_to_str(file_part[0])}_{flt_to_str(file_part[1])}_{flt_to_str(file_part[2])}_{flt_to_str(file_part[3])}_{flt_to_str(file_part[4])}_"
    file_parts_display_8_mp4.append(f"{temp_str}.mp4")
    file_parts_display_8_png.append(f"{temp_str}.png")

#print(file_parts_display_8_mp4, file_parts_display_8_png, file_parts_display_9_png, file_parts_display_9_mp4, sep=f'\n{"="*10}\n')

#print(f"file_parts_display_1_png of length {len(file_parts_display_1_png)}:", *file_parts_display_1_png, sep = '\n')
#exit()

#print("file_parts_display_5_mp4 =", file_parts_display_5_mp4)
titles = open("result_times.txt").readlines()
display_1_titles = ["Base Case"] + [i.strip('\n') for i in titles[:10]]
display_2_titles = ["Base Case"] + [i.strip('\n') for i in titles[10:20]]
display_3_titles = ["Base Case"] + [i.strip('\n') for i in titles[20:30]]
display_4_titles = ["Base Case"] + [i.strip('\n') for i in titles[30:40]]
display_5_titles = ["Base Case, A = 1.4"] + [i.strip('\n') for i in titles[40:49]]
display_6_titles = ["Base Case"] + [i.strip('\n') for i in titles[49:89]]
display_7_titles = ["Base Case"] + [i.strip('\n').replace(', =', ' =') for i in titles[89:99]]
display_8_titles = ["Base Case, A = 1.1"] + [i.strip('\n').replace(', =', ' =') for i in titles[100:]]
display_9_titles = [display_8_titles[-1], r"\(U_{\text{patched}}(x,\, A = 1,\, \mathcal{A}_0 = 1,\, \Omega = 0.2; \, \rho = 0.2)\) from (A = 1.080, b = 0.800, Ω = 0.230) \(\,\rightarrow \,\mathcal{L} = 1.4770857 \times 10^{-2}\)"]

print(*display_1_titles, sep='\n', end="\n\n")
print(*display_2_titles, sep='\n', end="\n\n")
print(*display_3_titles, sep='\n', end="\n\n")
print(*display_4_titles, sep='\n', end="\n\n")
print(*list(zip(range(1,1+len(display_5_titles)), display_5_titles)), sep='\n', end="\n\n")
print(*display_6_titles, sep='\n', end="\n\n")
print(*display_7_titles, sep='\n', end="\n\n")
print(*display_8_titles, sep='\n', end="\n\n")
print(*display_9_titles, sep='\n', end="\n\n")

display_1_title_tokens = [i.split() for i in display_1_titles[1:]]
display_2_title_tokens = [i.split() for i in display_2_titles[1:]]
display_3_title_tokens = [i.split() for i in display_3_titles[1:]]
display_4_title_tokens = [i.split() for i in display_4_titles[1:]]
display_5_title_tokens = [i.split() for i in display_5_titles[1:]]
display_6_title_tokens = [i.split() for i in display_6_titles[1:]]
display_7_title_tokens = [i.split() for i in display_7_titles[1:]]
display_8_title_tokens = [i.split() for i in display_8_titles[1:]]

display_1_titles = ["Base Case"]
display_2_titles = ["Base Case"]
display_3_titles = ["Base Case"]
display_4_titles = ["Base Case"]
display_5_titles = ["Base Case, A = 1.4"]
display_6_titles = ["Base Case"]
display_7_titles = ["Base Case"]
display_8_titles = ["Base Case, A = 1.1"]

print("display_1_title_tokens\n======================",*display_1_title_tokens, sep='\n', end="\n\n")
print("display_2_title_tokens\n======================",*display_2_title_tokens, sep='\n', end="\n\n")
print("display_3_title_tokens\n======================",*display_3_title_tokens, sep='\n', end="\n\n")
print("display_4_title_tokens\n======================",*display_4_title_tokens, sep='\n', end="\n\n")
print("display_5_title_tokens\n======================",*display_5_title_tokens, sep='\n', end="\n\n")
print("display_6_title_tokens\n======================",*display_6_title_tokens, sep='\n', end="\n\n\n")
print("display_7_title_tokens\n======================",*display_7_title_tokens, sep='\n', end="\n\n\n")
print("display_8_title_tokens\n======================",*display_8_title_tokens, sep='\n', end="\n\n\n")

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

for token_title in display_5_title_tokens:
    temp_str = f"Time spent learning {token_title[2]} = {token_title[6]} from {token_title[2]} = {token_title[4]} = {token_title[8]} seconds"
    if "not" in token_title:
        temp_str += r", \(\mathcal{L} = " f"{float(token_title[-1]):.4f}" r"\)"
    else:
        temp_str += r", \(\mathcal{L} < 1.3 \times 10^{-2}\)"
    display_5_titles.append(temp_str)

for token_title in display_6_title_tokens:
    token_title[6] = float(token_title[6])
    token_title[11] = float(token_title[11])
    token_title[4] = float(token_title[4])
    token_title[9] = float(token_title[9])
    temp_str = f"Time spent learning ({token_title[2]} = {token_title[6]:.3f}, {token_title[7]} = {token_title[11]:.3f}) from ({token_title[2]} = {token_title[4]:.3f}, {token_title[7]} = {token_title[9]:.3f}) = {token_title[13]} seconds"
    if "not" in token_title:
        temp_str += r", \(\mathcal{L} = " f"{float(token_title[-1]):.4f}" r"\)"
    else:
        temp_str += r", \(\mathcal{L} < 1.3 \times 10^{-2}\)"
    display_6_titles.append(temp_str)

for token_title in display_7_title_tokens:
    token_title[6] = float(token_title[6].strip(','))
    token_title[11] = float(token_title[11].strip(','))
    token_title[4] = float(token_title[4])
    token_title[9] = float(token_title[9])
    temp_str = f"Time spent learning ({token_title[2]} = {token_title[6]:.3f}, {token_title[7]} = {token_title[11]:.3f}) from ({token_title[2]} = {token_title[4]:.3f}, {token_title[7]} = {token_title[9]:.3f}) = {token_title[13]} seconds"
    if "not" in token_title:
        temp_str += r", \(\mathcal{L} = " f"{float(token_title[-1]):.4f}" r"\)"
    else:
        temp_str += r", \(\mathcal{L} < 1.3 \times 10^{-2}\)"
    display_7_titles.append(temp_str)

for token_title in display_8_title_tokens:
    token_title[6] = float(token_title[6].strip(','))
    token_title[11] = float(token_title[11].strip(','))
    token_title[4] = float(token_title[4])
    token_title[9] = float(token_title[9])
    token_title[14] = float(token_title[14])
    token_title[16] = float(token_title[16])
    token_title[18] = float(token_title[18])
    temp_str = f"Time spent learning ({token_title[2]} = {token_title[6]:.3f}, {token_title[7]} = {token_title[11]:.3f}, {token_title[12]} = {token_title[16]:.3f}) from ({token_title[2]} = {token_title[4]:.3f}, {token_title[7]} = {token_title[9]:.3f}, {token_title[12]} = {token_title[14]:.3f}) = {token_title[18]} seconds"
    if "not" in token_title:
        temp_str += r", \(\mathcal{L} = " f"{float(token_title[-1]):.4f}" r"\)"
    else:
        temp_str += r", \(\mathcal{L} < 1.3 \times 10^{-2}\)"
    display_8_titles.append(temp_str)

print(*display_1_titles, sep='\n', end='\n\n')
print(*display_2_titles, sep='\n', end='\n\n')
print(*display_3_titles, sep='\n', end='\n\n')
print(*display_4_titles, sep='\n', end='\n\n')
print(*display_5_titles, sep='\n', end='\n\n')
print(*display_6_titles, sep='\n', end='\n\n')
print(*display_7_titles, sep='\n', end='\n\n')
print(*display_8_titles, sep='\n', end='\n\n')
display_9_titles[0] = display_8_titles[-1]
print(*display_9_titles, sep='\n', end='\n\n')
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
        .scroll-container {
            overflow-x: auto;
            white-space: nowrap;
            width: 100%;
            max-width: 100%;
        }
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
<div class = "scroll-container" style = "text-align: center">
    <div class="latex-equation">
        $$\begin{align*}
        m \ddot{x} &= -\frac{d}{dx} \left(V_{\text{MT}}(x) + V_{\text{sech}}(x-\xi) \right) \\
        V_{\text{MT}}(x) &= \frac{1}{2} \Omega^2 x^2 \\
        V_{\text{sech}}(x-\xi) &= A \operatorname{sech}^{2}\left(b \left(x - {\xi}\right)\right)
        \end{align*}$$<br>
        $$
            U_{\text{patched}}(x,\, A,\, \mathcal{A}_0,\, \Omega; \rho) = \left\{
            \begin{array}{cc}
                \text{$U_{\text{eff}}(x,\, A,\, \mathcal{A}_0,\, \Omega; \rho)\, = - \dfrac{256 A \mathcal{A}_0^{2} x}{\left(e^{2 \mathcal{A}_0 x} - 1\right)^{5}} - \dfrac{256 A \mathcal{A}_0 \left(2 \mathcal{A}_0 x - 1\right)}{\left(e^{2 \mathcal{A}_0 x} - 1\right)^{3}} - \dfrac{128 A \mathcal{A}_0 \left(5 \mathcal{A}_0 x - 1\right)}{\left(e^{2 \mathcal{A}_0 x} - 1\right)^{4}} - \dfrac{32 A \mathcal{A}_0 \left(12 \mathcal{A}_0 x - 13\right)}{3 \left(e^{2 \mathcal{A}_0 x} - 1\right)^{2}} + \dfrac{32 A \mathcal{A}_0}{3 \left(e^{2 \mathcal{A}_0 x} - 1\right)} + \dfrac{2 \Omega^{2} \mathcal{A}_0 x^{2}}{3} $} & \quad\text{if } |\mathcal{X} | \geq \rho, \\[0.2cm]
                \text{$U_{\text{Taylor}}(x,\, A,\, \mathcal{A}_0,\, \Omega; \rho)\, = - \dfrac{64 A \mathcal{A}_0^{3} x^{2}}{105} + \dfrac{16 A \mathcal{A}_0}{15} + \dfrac{2 \Omega^{2} \mathcal{A}_0 x^{2}}{3}$} & \quad\text{if } |\mathcal{X} | < \rho.
            \end{array}\right.
        $$
    </div>
</div>
<h3><u>Constraints</u></h3>
<div class = "scroll-container" style = "text-align: center">

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
</div>
<h3 class="latex-equation underline">Loss Function \(\mathcal{L}\)</h3>
<div class = "scroll-container" style = "text-align: center">

    <div class="latex-equation underline">
        $$\begin{align*}
        &\overrightarrow{\text{MSE}} \equiv \epsilon\cdot \left( \left(x^* - \vec{x}\right)\odot\left(x^* - \vec{x}\right)\right) + \zeta\cdot\left(\vec{v}\odot\vec{v}\right) + \eta\left(\left(x^* - \vec{\xi}\right)\odot\left(x^* - \vec{\xi}\right)\right) \\
        &\text{MSE}_{\text{best}} = \min\left(\overrightarrow{\text{MSE}}\right) \\
        &t_{\text{best}} = \underset{t_i\, \in \,\vec{t}}{\text{argmin}} \left(\overrightarrow{\text{MSE}}\right) \\
        &\text{Smoothness Penalty} = \frac{1}{t_{\text{best}}}\sum_{t=1}^{t_{\text{best}}} \left(\xi_{t} - \xi_{t-1}\right)^2 \\
        &\text{Max}_{v} = \max\left(\left|\vec{v}_{[t=1:t_{\text{best}}]}\right|\right) \\
        &\text{Max}_{\xi} = \max\left(\left|\vec{\xi}_{[t=1:t_{\text{best}}]}\right|\right) \\[0.5cm]
        &\boxed{\mathcal{L} = \text{MSE}_{\text{best}} + \alpha \cdot \text{Smoothness Penalty} + \beta \cdot t_{\text{best}} + \gamma \cdot \text{Max}_{v} + \delta \cdot \text{Max}_{\xi}}
        \end{align*}$$
        <br> Where:
        $$\begin{align*}
        &\epsilon = \zeta = \eta = 1 \\
        &\alpha = \beta = \gamma = \delta = 10^{-3} \\
        &\vec{x} = [x_1, x_2, ..., x_{T_{\text{Th}}}] \\
        &\vec{v} = [v_1, v_2, ..., v_{T_{\text{Th}}}] \\
        &\vec{\xi} = [\xi_1, \xi_2, ..., \xi_{T_{\text{Th}}}] 
        \end{align*}$$
    </div>
</div>
<br>
<h3><u>Displays</u></h3>
"""
    for i, display in enumerate(file_data):
        first_png = display["png"][display["index"]]
        first_mp4 = display["mp4"][display["index"]]
        html_content += f"""
        <h4 id="title{i}" class="title latex-equation scroll-container">{file_data[i]['titles'][file_data[i]['index']]}</h4> 
        <div class="display" id="display{i}">
            <button onclick="cycleDisplay({i}, -1)">&#8592;</button>
            <img id="img{i}" src="{first_png}" alt="Display {i + 1} Image">
            <video id="video{i}" src="{first_mp4}" controls></video>
            <button onclick="cycleDisplay({i}, 1)">&#8594;</button>
        </div>
        """

    html_content += r"""
    <h4 id="A_1_point_5_from_scratch" class="title latex-equation scroll-container">\(\mathbf{A = 1.5\,}\) from scratch \(\mathbf{\rightarrow \mathcal{L} = 0.0144}\)</h4> 
    <div class="display">
        <img id="img_A_1_point_5_from_scratch" src="trajectory_data_IC_2_point_3981_1_point_5_1_point_0_1_point_0_0_point_2_random_initialized.png" alt="A = 1.5 from Scratch">
        <video id="movie_A_1_point_5_from_scratch" src="trajectory_data_IC_2_point_3981_1_point_5_1_point_0_1_point_0_0_point_2_random_initialized.mp4" controls></video>
    </div>"""
    
    html_content += r"""
    <h4 id="A_point_0675_B_point_666_from_A_1_B_1" class="title latex-equation scroll-container">\(\mathbf{A = 0.0675, \:b = 0.666}\) from \(\mathbf{A = 1, \:b = 1}\) \(\mathbf{\rightarrow \mathcal{L} < 0.013}\)</h4> 
    <div class="display">
        <img id="img_A_point_0675_B_point_666_from_A_1_B_1" src="trajectory_data_IC_0_point_848359_0_point_067475_0_point_665792_1_point_0_0_point_2_.png" alt="A = 0.0675, B = 0.666 from Scratch">
        <video id="movie_A_point_0675_B_point_666_from_A_1_B_1" src="trajectory_data_IC_0_point_848359_0_point_067475_0_point_665792_1_point_0_0_point_2_.mp4" controls></video>
    </div>"""
    
    html_content += r"""
    <h4 id="Paul_From_A_point_0675_B_point_666" class="title latex-equation scroll-container">\(\mathbf{V_{\text{eff}}(x, A = 1, B = 0.1) = \frac{8B\left(Ae^{4A(x-\xi)}(x-\xi)+Ae^{2A(x-\xi)}(x-\xi)-e^{4A(x-\xi)}+e^{2A(x-\xi)}\right)}{e^{6A(x-\xi)}-3e^{4A(x-\xi)}+3e^{2A(x-\xi)}-1}}\:\) from \(\mathbf{A = 0.0675, \:b = 0.666}\) \(\mathbf{\rightarrow \mathcal{L} = 0.016351}\)</h4> 
    <div class="display">
        <img id="img_Paul_From_A_point_0675_B_point_666" src="trajectory_data_IC_0_point_787127_1_0_point_1_1_point_0_0_point_2_Paul_.png" alt="A = 0.0675, B = 0.666 from Scratch">
        <video id="movie_Paul_From_A_point_0675_B_point_666" src="trajectory_data_IC_0_point_787127_1_0_point_1_1_point_0_0_point_2_Paul_.mp4" controls></video>
    </div>"""

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
    {"png": file_parts_display_5_png, "mp4": file_parts_display_5_mp4, "index": 0, "titles": display_5_titles},
    {"png": file_parts_display_6_png, "mp4": file_parts_display_6_mp4, "index": 0, "titles": display_6_titles},
    {"png": file_parts_display_7_png, "mp4": file_parts_display_7_mp4, "index": 0, "titles": display_7_titles},
    {"png": file_parts_display_8_png, "mp4": file_parts_display_8_mp4, "index": 0, "titles": display_8_titles},
    {"png": file_parts_display_9_png, "mp4": file_parts_display_9_mp4, "index": 0, "titles": display_9_titles}
]


#[os.system(f"ls {i}") for i in file_parts_display_8_png]
#exit()

# Generate and save HTML
html_content = generate_html(file_data)
with open("control_displays.html", "w") as f:
    f.write(html_content)

print("HTML file generated as 'control_displays.html'")
