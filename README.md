# KleoProj
Personal project to recreate the periodic orbits from [A. Abad et al., 2024](https://doi.org/10.1016/j.asr.2024.10.017)

This project uses KleoPy: a custom python package implemented to help find periodic orbits around *216-Kleopatra Asteroid*.

## Kleopy
Custom tools for 216-Kleopatra using the Dipole-Segment model.

*Submodules:*
1. constants : Physical constants and parameters for 216-Kleopatra
2. grid_search : Numerical methods for finding symmetric periodic orbits
3. integrator : Numerical integrators for solving ODEs
4. orbital_eq : Equations related to the orbital dynamics of 216-Kleopatra

## License
[BSD 3-Clause License](LICENSE)

## Authors
- Gavin Athaya S. [(@gavinathaya)](https://github.com/gavinathaya)
- Fernanda [(@fernandaputra21)](https://github.com/fernandaputra21)

## Acknowledgements
*This whole project and kleopy leverages the following python packages*
1. [NumPy](https://numpy.org)
2. [SciPy](https://scipy.org)
3. [Matplotlib](https://matplotlib.org)

## Bibliography
1. [Alberto Abad, Antonio Elipe, Alessandra F.S. Ferreira, Periodic orbits around 216-Kleopatra asteroid modelled by a dipole-segment, *Advances in Space Research*, Volume 74, Issue 11, 2024, Pages 5687-5697, ISSN 0273-1177.](https://doi.org/10.1016/j.asr.2024.10.017)
2. [Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming with NumPy. *Nature* 585, 357–362 (2020). DOI: 10.1038/s41586-020-2649-2.](https://doi.org/10.1038/s41586-020-2649-2)
3. [Pauli Virtanen, Ralf Gommers, Travis E. Oliphant, Matt Haberland, Tyler Reddy, David Cournapeau, Evgeni Burovski, Pearu Peterson, Warren Weckesser, Jonathan Bright, Stéfan J. van der Walt, Matthew Brett, Joshua Wilson, K. Jarrod Millman, Nikolay Mayorov, Andrew R. J. Nelson, Eric Jones, Robert Kern, Eric Larson, CJ Carey, İlhan Polat, Yu Feng, Eric W. Moore, Jake VanderPlas, Denis Laxalde, Josef Perktold, Robert Cimrman, Ian Henriksen, E.A. Quintero, Charles R Harris, Anne M. Archibald, Antônio H. Ribeiro, Fabian Pedregosa, Paul van Mulbregt, and SciPy 1.0 Contributors. (2020) SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. *Nature Methods*, 17(3), 261-272. DOI: 10.1038/s41592-019-0686-2.](https://doi.org/10.1038/s41592-019-0686-2)
4. [J. D. Hunter, "Matplotlib: A 2D Graphics Environment", Computing in Science & Engineering, vol. 9, no. 3, pp. 90-95, 2007.](https://doi.org/10.5281/zenodo.14940554)