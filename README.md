# Main codes: active_polymer_adsorption.f

Langevin dynamics simulation of end-grafted active polymer chain. It was written in Fortran 90. Detailed simulation methods can be found in the following cites:

(1) M. B. Luo and Y. F. Shen, Langevin dynamics simulations for the critical adsorption of end-grafted active polymers, Soft Matter, 2024, 20, 5113–5121.

(2) Y. F. Shen and M. B. Luo, Knotting and adsorption of end-grafted active polymers, Soft Matter, submitted.

[github repo](https://github.com/mbluo-zju/active-polymer)


## Usage：
Windows:  active_polymer_adsorption.exe 

Linux:     ./active_polymer_adsorption.X 

## Format of input file "parameters_input.txt"

```
'n65k0f2Dr3s3'	output file names
0    'x.cfg'    if the first number is 1, then read in initial conformation

200		       Number of samples
60.0d0	60.0d0	40.0d0	system size x y z

1	65	number of polymer chain and length (the first bead is active bead)  
1.0d0	bead size
30.0d0		FENE strength
1.5d0	1.2d0	maximum FENE bond for passive part and the first bond connected the active bead
0.0d0   0.0d0   Bending strength and equilibrium bond angle

1.0d0	1.12d0 	WCA potential   

2.0d0		Active force Fa
3.0d0	1.0d0	Diffusion coefficients D_r and D_t
0.1d0      Rotational inertia J          

1.0d0	 temperature of solvent

0.005d0		timestep for Langevin dynamics 
2000	20000 	500	  500  Equilibrium time, Statistical time, Time for recording conformation, Time for output conformation

3         suface type
2.5d0      rcut for polymer-surface
1  36      Number of polymer-surface interaction Eps (starting No and ending No)    
0.05  0.1  0.15  0.2  0.3  0.4  0.5   0.6  0.65  0.7 0.75  0.8 0.85  0.9  0.95  1.0  1.05  1.1 1.15 1.2  1.25 1.3 1.35  1.4  1.45  1.5  1.6  1.7   1.8  1.9   2.0  2.2  2.4  2.6  2.8  3.0     Interaction strength for simulation
```  

Every energy has its own conformation file with suffix 00xx

## Output:

### Data file:

Format:
```
"Eps"  "mean-square end-to-end distance R^2"  "surface normal component of mean-square end-to-end distance R^2_z"  "mean-square radius of gyration R_G^2"  "surface normal component of mean-square radius of gyration R_G^2_z"  "mean height of mass center of polymer chain z_c"  "mean bond length"  "two-dimensional shape parameter A_2D"  "mean number of adsorbed monomer M"   "mean square M"   "mean interaction energy between polymer and surface E"   "mean square E"
```
 
### Conformation file:

Record coordinates of the whole polymer chain

Format:

```
serial number   polymer ID   monomer ID   x  y  z   monomer species
```
