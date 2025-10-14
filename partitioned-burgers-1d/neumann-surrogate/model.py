import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn.utils import weight_norm

def pad_with_ghost_cells(input_seq, bc_left, bc_right):
    return torch.cat([bc_left, input_seq, bc_right], dim=1)

class LinearExtrapolationPadding1D(nn.Module):
    """Applies 'same' padding using linear extrapolation."""
    def __init__(self, kernel_size: int, dilation: int = 1):
        super().__init__()
        self.pad_total = dilation * (kernel_size - 1)
        self.pad_beg = self.pad_total // 2
        self.pad_end = self.pad_total - self.pad_beg

    def forward(self, x):
        # Don't pad if not necessary
        if self.pad_total == 0:
            return x

        ghost_cell_left = x[:, :, :1]
        ghost_cell_right = x[:, :, -1:]

        # Calculate the gradient at each boundary
        grad_left = x[:, :, 1:2] - ghost_cell_left
        grad_right = ghost_cell_right - x[:, :, -2:-1]

        # Higher order finite difference gradient approximation
        # grad_left = ( -11 * ghost_cell_left + 18 * x[:, :, 1:2] - 9 * x[:, :, 2:3] + 2 * x[:, :, 3:4]) / 6
        # grad_right = (11 * ghost_cell_right - 18 * x[:, :, -2:-1] + 9 * x[:, :, -3:-2] - 2 * x[:, :, -4:-3]) / 6

        # 3. Extrapolated padding tensors
        left_ramp = torch.arange(self.pad_beg, 0, -1, device=x.device, dtype=x.dtype).view(1, 1, -1)
        left_padding = ghost_cell_left - left_ramp * grad_left

        right_ramp = torch.arange(1, self.pad_end + 1, device=x.device, dtype=x.dtype).view(1, 1, -1)
        right_padding = ghost_cell_right + right_ramp * grad_right
        
        return torch.cat([left_padding, x, right_padding], dim=2)

class ResidualBlock1D(nn.Module):
    """A residual block that uses custom 'same' padding with linear extrapolation and weight normalization."""
    def __init__(self, channels, kernel_size=3, activation=nn.ReLU):
        super(ResidualBlock1D, self).__init__()
        self.activation = activation()
        # Apply weight normalization 
        self.conv1 = weight_norm(nn.Conv1d(channels, channels, kernel_size, padding='valid', bias=True))
        self.ghost_padding1 = LinearExtrapolationPadding1D(kernel_size)
        self.conv2 = weight_norm(nn.Conv1d(channels, channels, kernel_size, padding='valid', bias=True))
        self.ghost_padding2 = LinearExtrapolationPadding1D(kernel_size)

    def forward(self, x):
        identity = x
        
        out = self.ghost_padding1(x)
        out = self.conv1(out)
        out = self.activation(out)
        
        out = self.ghost_padding2(out)
        out = self.conv2(out)

        return self.activation(out) + identity
    
class CNN_RES(nn.Module):
    """
    A CNN with residual blocks for 1D data.
    Expects a pre-padded input with ghost_cells//2 number ghost cells on each side.
    Applies a custom linear extrapolation padding for inner layers.
    """
    def __init__(self, hidden_channels, num_blocks=2, kernel_size=3, activation=nn.ReLU, ghost_cells=2):
        super(CNN_RES, self).__init__()
        self.activation = activation()
        self.hidden_channels = hidden_channels
        self.num_blocks = num_blocks
        self.kernel_size = kernel_size
        assert ghost_cells % 2 == 0, "ghost_cells must be even"
        self.ghost_cells = ghost_cells
        
        self.ghost_padding = LinearExtrapolationPadding1D(self.ghost_cells + self.kernel_size)

        # Apply weight normalization to the input convolution
        self.conv_in = weight_norm(nn.Conv1d(1, hidden_channels, kernel_size=1, bias=True))
        
        layers = [ResidualBlock1D(hidden_channels, kernel_size, activation=activation) for _ in range(num_blocks)]
        self.res_blocks = nn.Sequential(*layers)
        
        self.conv_out = nn.Conv1d(hidden_channels, 1, kernel_size=1)
    
    def forward(self, x):

        if x.dim() == 2:
            x = x.unsqueeze(1) # Add channel dim: (B, 1, L)

        if not self.ghost_cells == 0:
            x_padded = self.ghost_padding(x)
            
        else:
            x_padded = x

        total_pad_each_side = self.ghost_padding.pad_beg + self.ghost_cells // 2

        out = self.activation(self.conv_in(x_padded)) # no extra padding here
        out = self.res_blocks(out)
        out = self.conv_out(out) # no extra padding here

        if not self.ghost_cells == 0:
            out = out[:, :, total_pad_each_side:-total_pad_each_side] # remove ghost cells, return only internal domain
        
        return out.squeeze(1)