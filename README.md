## Ordering specification for POCRM with Triple Combinations
This repository contains code used for the paper _Applying the partial order continual reassessment method to high-dimensional treatment combinations_ by Weishi Chen, Li Liu, Nolan A. Wages, and Pavel Mozgunov.

## Structure
- The folder `R functions` contains R function files:
  - `SimScen.R` contains code to simulate scenarios under 3-dimensional grids.
  - `ListOrder.R` contains code to list out all possible orderings for the 12 combination example.
- The folder `RData` contains RDatas 
  - `148Orderings.RData` contains all 148 possible orderings for the 12 combinations example.
  - `OrderScen_12combo.RData` contains the 47 order-scenarios for the 12 combinations example.

## Dependencies
This project is written in R. Key R packages used (as seen in the code) include:
- dfcrm
- pocrm
- combinat
- tidyverse

##  License
MIT License. See `LICENSE` file for details.
