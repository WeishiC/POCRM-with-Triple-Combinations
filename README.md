## Ordering specification for POCRM with Triple Combinations
This repository contains code used for the paper _Applying the partial order continual reassessment method to high-dimensional treatment combinations_ by Weishi Chen, Li Liu, Nolan A. Wages, and Pavel Mozgunov.

## Structure
- The folder `R functions` contains R function files:
  - `SimScen.R` contains code to simulate scenarios under 3-dimensional grids.
  - `ListOrder.R` contains code to list out all possible orderings for the 12 combination example.
  - `Adding.R` contains the code to implement the *Adding* step of the *Adding-Refining algorithm*.
  - `POCRM.R` contains the code to implement and simulate under 2-stage likelihood-based POCRM (Wages et al, 2011).
- The folder `RData` contains RDatas 
  - `148Orderings.RData` contains all 148 possible orderings for the 12 combinations example.
  - `ScenAgnosticOrd.RData` contains 8 scenario-agnostic orderings for the 12 combinations example.
  - `ScenSpecOrd.RData` contains 4 scenario-agnostic orderings for the 12 combinations example.
  - `AROrd.RData` contains a list of 10 sets of orderings obtained from the Adding-Refining algorithm for the 12 combinations example.
  - `InconsisOrd.RData` contains an inconsistent choice of 4 orderings for the 12 combinations example.
  - `OrderScen_12combo.RData` contains the 47 order-scenarios for the 12 combinations example.
  - `SimScen.RData` contains the 12 scenarios used for the simulation study in the paper.

## Dependencies
This project is written in R. Key R packages used (as seen in the code) include:
- dfcrm
- pocrm
- nnet
- combinat
- tidyverse

## References
- Chen, W., Liu, L., Wages, N.A. and Mozgunov, P. *Applying the partial order continual reassessment method to high-dimensional treatment combinations.* 2025
- Wages, N.A., Conaway, M.R. and O'Quigley, J. *Dose-finding design for multi-drug combinations.* Clinical Trials 2011; 8: 380-389.

##  License
MIT License. See `LICENSE` file for details.
