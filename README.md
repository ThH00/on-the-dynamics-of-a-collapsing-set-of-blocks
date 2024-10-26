# stacked-blocks-2d
This code simulates the motion of a vertical stack of two dimensional blocks.
It accompanies the paper ....TBD....

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


## Debugging and code amelioration

### Value of damping factor $\rho_{\infty}$
- Date: 03/11/2024
- Problem: Some values of $\rho_{\infty}$ lead to a singular Jacobian. Refer to [Capobianco paper](https://doi.org/10.1002/nme.6801) for more information about $\rho_{\infty}$.
- Solution: I implement a 'try: except:' loop that goes though $\rho_{\infty}$ values until the code converges. Changing the damping factor values should not affect the accuracy of the results because I am using at iterative method at each time step and the code moves on to the next time step only when the desired error tolerance is achieved.
- Future work: Explore the details of the effect of $\rho_{\infty}$ on the singularity of the Jacobian.
- [x] Resolved

### Infinite loop of incrementation of damping factor
- Date: 04/01/2024
- Problem: To solve the previous problem, we resort to incrementing $\rho_{\infty}$ which opens up the possibility of an infinite loop.
- Solution: If `TimeoutError`, fix contact regions and check for convergence. If fixing contact region does not work, increment $\rho_{\infty}$.
- Future work: Maybe I can save `nu` iteration error and make corrolation between cyclical error and increasing error and the course of action to be taken.
- [x] Resolved

### Changing contact regions within time step
- Date: 03/11/2024
- Problem: Code does not converge because it iterates between two or more contact regions configurations.
- Solution: There might be bifurcations here. Fix a contact region and explore code convergence.
- [x] Resolved on 04/12/2024

### Error tolorance
- Date: 03/15/2024
- Problem: Errors to the order e-6 in acceleration compound in q and u and lead to impact (sudden change in contact regions) that should not happen.
- Solution: Decrease error tolerance to e-10. I also though that nondimensionalization should help with this.
- Future work: Concern of computer accuracy arises.
- Comments: This problem might have been due to a bug that is now fixed.
- [x] Resolved

### Nondimensionalization
- Date: 03/15/2024
- Problem: Descripancies in the values of the numerical entries in the Jacobian. Some values of length parameter lead to to a poorly scalted Jacobian (not invertible) because determinant is not zero, but is too small (close to machine epsilon).
- Solution: Rescale to make all terms in the Jacobian of the same order.
- [ ] Resolved

## Elaboration on the implementation of bifurcation exploration

### `def solve`
- create a new function: `solve`
- the first time `solve` is run, no contact region is given to the code
- any time `solve` is ran after that, it starts from a bifurcation point
- to save the results, (1) keep track of the parent of each bifurcation (see `num_bif_contacts`) (2) save the full displacement array for each bifurcation.

### `num_bif_contacts`
- is an array with dimensionas number_of_leaves x 4
- the first column save the leave number at which the bifurcation happens (`leaves_counter`)
- the second column saves the iteration at which bifurcation happens (`iter_bif`)
- the third column saves the number of unique contact region per bifuraction (`n_unique_contacts`)
- the fourth column saves the converged number of bifurcations out of the total number of unque contacts (that number is initialized to `n_unique_contacts` and is then decremented by 1 everytime I arrive to the else statement is `def solve_bifurcation`).


