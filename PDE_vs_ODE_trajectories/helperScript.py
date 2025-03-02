from glob import glob
from os import system
delimiter = ", "
mass_images = [i.split('_') for i in glob("mass*png")]
trajectory_images = [i.split('_') for i in glob("trajectory*png")]
mass_images.sort(key = lambda x: int(x[-1][:x[-1].find('.')]))
trajectory_images.sort(key = lambda x: int(x[-1][:x[-1].find('.')]))
mass_images = ['"'+'_'.join(i)+'"' for i in mass_images]
trajectory_images = ['"'+'_'.join(i)+'"' for i in trajectory_images]
titles = [rf'"\\(v_0 = {i}\\)"' for i in range(1, 11)]
#[system(f"ls {i}") for i in mass_images]
#[system(f"ls {i}") for i in trajectory_images]
print(*mass_images, sep=delimiter)
print("="*40)
print(*trajectory_images, sep=delimiter)
print(*titles, sep=delimiter)
