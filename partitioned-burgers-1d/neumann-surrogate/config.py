import torch

# Model architecture
INPUT_SIZE = 128 + 2 # +2 for ghost cells
HIDDEN_SIZE = 64 # num filters
OUTPUT_SIZE = 128

assert INPUT_SIZE >= OUTPUT_SIZE, "Input size must be greater or equal to output size."
assert (INPUT_SIZE - OUTPUT_SIZE) % 2 == 0, "Input and output sizes must differ by an even number (for ghost cells)."

NUM_RES_BLOCKS = 4
KERNEL_SIZE = 5
ACTIVATION = torch.nn.ReLU

MODEL_NAME = "CNN_RES_UNROLL_7.pth"