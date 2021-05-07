# BDEtools

An open-source MATLAB package for solving sets of Boolean delay equations

## Introduction

Systems of Boolean Delay Equations (BDEs) – in which time is continuous but state is binary – are capable of generating surprisingly complex behaviour, despite their apparent simplicity [1,2]. In addition to simulating convergence to steady states, BDEs can also generate periodic and quasiperiodic oscillations, 𝑚:𝑛 frequency locking and even chaotic dynamics [3]. Furthermore, the enumerability of the Boolean update functions and the compact parametrisation resulting from discretisation means that BDE systems can be readily leveraged to generate low-level descriptions of physical systems, from which more quantitative model formulations (e.g. differential equations) can be constructed [4]. The utility of BDE modelling in this regard has been demonstrated in several fields, including computational biology and climate science, but the use of BDEs is still primarily restricted to a few research laboratories. In order to facilitate the wider adoption of the BDE formalism by the computational modelling community, BDEtools has been developed to enable researchers to construct, solve and analyse BDE systems in a straightforward fashion.

## Overview of file structure

<ul> 
<li> code/ -- BDEtools solvers and utility routines
<li> examples/ -- three examples illustrating the use of the solvers and utility routines: (i) a simple two-component model that generates accumulating switches; (ii) a climate model of the El-Ni&ntildeo Southern Oscillation; and (iii) a computational biology model of circadian oscillations in the fungus <i>N. crassa</i>  
<li> models/ -- BDEtools implementations of the circadian clock models presented in [1].  
</ul>  



## References

[1] D. Dee and M. Ghil. Boolean difference equations, i: Formulation and dynamic behavior. SIAM J. Appl. Math., 44(1):111–126, 1984.
[2] M. Ghil and A. Mullhaupt. Boolean delay equations. ii. periodic and ape- riodic solutions. J. Stat. Phys., 41(1):125–173, 1985.
[3] M. Ghil, I. Zaliapin, and B. Coluzzi. Boolean delay equations: A simple way of looking at complex systems. Physica D, 237(23):2967–2986, 2008.
[4] O.E. Akman, S. Watterson, A. Parton, N. Binns, A.J. Millar, and P. Ghazal. Digital clocks: simple Boolean models can quantitatively de- scribe circadian systems. J. Roy. Soc. Interface, 9(74):2365–2382, 2012.
