# BDEtools

An open-source MATLAB package for solving sets of Boolean delay equations

## Introduction

Systems of Boolean Delay Equations (BDEs) – in which time is continuous but state is binary – are capable of generating surprisingly complex behaviour, despite their apparent simplicity. In addition to simulating convergence to steady states, BDEs can also generate periodic and quasiperiodic oscillations, 𝑚:𝑛 frequency locking and even chaotic dynamics. Furthermore, the enumerability of the Boolean update functions and the compact parametrisation resulting from discretisation means that BDE systems can be readily leveraged to generate low-level descriptions of physical systems, from which more quantitative model formulations (e.g. differential equations) can be constructed. The utility of BDE modelling in this regard has been demonstrated in several fields, including computational biology and climate science, but the use of BDEs is still primarily restricted to a few research laboratories. In order to facilitate the wider adoption of the BDE formalism by the computational modelling community, BDEtools has been developed to enable researchers to construct, solve and analyse BDE systems in a straightforward fashion.

## Overview of file structure

\ul
