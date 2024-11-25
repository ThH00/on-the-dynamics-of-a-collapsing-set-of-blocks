# stacked-blocks-2d
This code simulates the motion of a vertical stack of two dimensional blocks.
It accompanies the paper *On the Dynamics of a Collapsing Set of Blocks* submitted to publication in *Acta Mechanica*.

## Features of the code.
- The motion of the bottom block is prescribed.
- Forces and moments can be applied to any block.

## branch-by-branch approach vs iteration-by-iteration approach

To solve for the motion of a stack of blocks produced by the oscillation of the bottom block, two approaches were followed. 

For the first approach, which we refer to as branch-by-branch approach, an initial solution branch is solved while recording the points along the branch where other solutions are possible. The code then revisits these potential bifurcation points and computes the additional branches. This method first resolves solutions arising at the end of the simulation and then progresses backward to address potential bifurcations occurring earlier in the process.

In the second approach, which we refer to as iteration-by-iteration approach, we compute all possible branches at each step of the numerical scheme. For a given branch, if none of the possible solutions converge, even after fixing the contact states and adjusting the parameter $\rho_{\infty,0}$, the branch is discarded.
 

The branch-by-branch approach was used to compute the solutions shown in Figure 6, while the iteration-by-iteration approach was employed to identify the distinct solutions in Figure 5 that occurred  earlier in the simulations.

> [!NOTE]
> The blocks are numbered from the bottom up. The bottom block with prescribed motion is block 0. The first block above it is block 1, etc.

## Instructions.
- To reproduce the results from Section 4, run `main-lire-tower.py`
- To reproduce the results from Section 5, run `runner-oscillation-of-bottom-block.py` for Figure 6 and `main-oscillation-of-bottom-block.bbb` for Figure 5.
- To reproduce the results from Section 6, run `runner-sliding-bottom-block.py`

You may indicate your output folder in the corresponding Matlab `.m` scripts and visualize the results.