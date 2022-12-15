# Griffiths & Jenkins (2022)
Matlab code to recreate the simulation study in "An estimator for the recombination rate from a continuously observed diffusion of haplotype frequencies".

## Citation
Jenkins, P. A., and Griffiths, R. C. An estimator for the recombination rate from a continuously observed diffusion of haplotype frequencies. In preparation.

## License information
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

## Description of files

driftmu.m - Function to compute the drift mu_ij on p13 of the paper. Input: a 2x2 matrix x and (optional) arguments rho, theta_A, theta_B, P_A, P_B. Returns the recombination and mutation components of mu_ij separately.

diffusionsigma.m - Function to compute the diffusion sigma on p19 of the paper. Input: a 2x2 matrix x.

diffusionpathEuler.m - Function to implement the Euler-Maruyama method. Input: parameters theta_A, theta_B, P_A, P_B, rho, R (number of replicates), N (number of timesteps), dt (stepsize), x_0. Returns a table containing, for each replicate: rhohat,rhoerror,minx (minimum of X_{ij}(t) across (i,j,t),rhohatN (numerator of rho-hat),rhoerrorN (numerator of rho-error),I,Iinf (indicator for I_T = infinity).

main.m - Wrapper code to run simulation study in reported in the paper. Outputs Table 1 (up to random seed).

results_tA*_rho*_N1e6_x0b.txt - Output from main.m used for Table 1 in the paper.

table_tA*_N1e6_x0b.txt - Table 1 in the paper.

## Version history

v1.0 (2022-12-15) --- First release
