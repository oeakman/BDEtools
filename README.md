# BDEtools

An open-source MATLAB package for solving sets of Boolean delay equations

## Introduction

Systems of Boolean Delay Equations (BDEs) – in which time is continuous but state is binary – are capable of generating surprisingly complex behaviour, despite their apparent simplicity [1,2]. In addition to simulating convergence to steady states, BDEs can also generate periodic and quasiperiodic oscillations, 𝑚:𝑛 frequency locking and even chaotic dynamics [3]. Furthermore, the enumerability of the Boolean update functions and the compact parametrisation resulting from discretisation means that BDE systems can be readily leveraged to generate low-level descriptions of physical systems, from which more quantitative model formulations (e.g. differential equations) can be constructed [4]. The utility of BDE modelling in this regard has been demonstrated in several fields, including computational biology and climate science [4-10], but the use of BDEs is still primarily restricted to a few research laboratories. In order to facilitate the wider adoption of the BDE formalism by the computational modelling community, BDEtools has been developed to enable researchers to construct, solve and analyse BDE systems in a straightforward fashion.

## Overview of file structure

<ul> 
<li> <b>code/</b> -- BDEtools solvers and utility routines
<li> <b>examples/</b> -- three examples illustrating the use of the solvers and utility routines: 
<ol>  
  <li> a simple two-component model that generates accumulating switches [2]
  <li> a climate model of the El-Ni&ntildeo Southern Oscillation [6] 
  <li> a computational biology model of circadian oscillations in the fungus <i>N. crassa [4] </i>  
</ol>  
<li> <b>models/</b> -- BDEtools implementations of the circadian clock models presented in [4].
<li> <b>models/</b> -- unit tests to check the correct impementations of the systems in <b>models/</b>
</ul>  

## License

BDEtools is released under the MIT license.

## References

[1] D. Dee and M. Ghil. Boolean difference equations, i: Formulation and dynamic behavior. SIAM J. Appl. Math., 44(1):111–126, 1984. 
<br>
[2] M. Ghil and A. Mullhaupt. Boolean delay equations. ii. periodic and ape- riodic solutions. J. Stat. Phys., 41(1):125–173, 1985.
<br>
[3] M. Ghil, I. Zaliapin, and B. Coluzzi. Boolean delay equations: A simple way of looking at complex systems. Physica D, 237(23):2967–2986, 2008.
<br>
[4] O.E. Akman, S. Watterson, A. Parton, N. Binns, A.J. Millar, and P. Ghazal. Digital clocks: simple Boolean models can quantitatively de- scribe circadian systems. J. Roy. Soc. Interface, 9(74):2365–2382, 2012.
<br>
[5] M. Ghil, A. Mullhaupt, and P. Pestiaux. Deep water formation and quaternary glaciations. Climate Dynamics, 2(1):1–10, 1987.
<br>
[6] A. Saunders and M. Ghil. A Boolean delay equation model of ENSO variability. Physica D, 160:54–78, 2001.
<br>
[7] H. Öktem, R. Pearson, and K. Egiazarian. An adjustable aperiodic model class of genomic interactions using continuous time Boolean networks (Boolean delay equations). Chaos, 13(4):1167–1174, 2003.
<br>
[8] R. Zhang, H.L.D. de S.Cavalcante, Z. Gao, D.J. Gauthier, J.E.S. Socolar, M.M. Adams, and D.P. Lathrop. Boolean chaos. Phys. Rev. E, 80:045202, 2009.
<br>
[9] K. Doherty, K. Alyahya, O.E. Akman, and J.E. Fieldsend. Optimisation and landscape analysis of computational biology models: A case study. In Proc. GECCO 2017, pages 1644–1651, 2017.
<br>
[10] O.E. Akman and J.E. Fieldsend. Multi-objective optimisation of gene regulatory networks: Insights from a Boolean circadian clock model. In Proc. BICOB 2020, volume 70, pages 149–162, 2020.

<hr>
&#169; Akman Laboratory of Automated Biotechnology, 2021
