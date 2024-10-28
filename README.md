# stacked-blocks-2d
This code simulates the motion of a vertical stack of two dimensional blocks.
It accompanies the paper *On the Dynamics of a Collapsing Set of Blocks* submitted to publication in *Acta Mechanica*.

## Features of the code.
- The motion of the bottom block is prescribed.
- Forces and moments can be applied to any block.

> [!NOTE]
> The blocks are numbered from the bottom up. The bottom block with prescribed motion is block 0. The first block above it is block 1, etc.

## Instructions.
- To reproduce the results from Section 4, run `main-lire-tower.py`
- To reproduce the results from Section 5, run `runner-oscillation-of-bottom-block.py`
- To reproduce the results from Section 6, run `runner-sliding-bottom-block.py`

You may indicate your output folder in the corresponding Matlab `.m` scripts and visualize the results.