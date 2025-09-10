import numpy as np
from scipy.io import savemat, loadmat

input_arr = loadmat("output.mat")["out"]
output_arr = input_arr + 5
output_filepath = "output.mat"
savemat(output_filepath, {"out": output_arr})
print(f"file created {output_filepath}")
print(loadmat(output_filepath)["out"])
