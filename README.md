# Supplementary scripts for "Spectral Invariants of Circular Tensor Networks on the Moduli Space of Fano 3-Folds"

## Requirements
- Python 3.10+
- mpmath >= 1.3.0
- numpy >= 1.24.0
- scipy >= 1.10.0
- matplotlib >= 3.7.0

## Installation
pip install mpmath numpy scipy matplotlib

## Scripts and usage

- Alpha_MPS_ClosedForm_MonteCarlo.py : Closed-form evaluation of α⁻¹, MPS convergence (Table 1), Monte Carlo sampling (Table 2)
- PEPS_5D_NODE_SIMULATOR_Quantum_Gravity.py : VEV saturation, hinge integration, topological invariant -40 (Appendix C)
- Isomorphism_PF_DS_Fano_2-22.py : Minkowski period sequence, factorizations (Section 6)
- Algebraic_Collapse_Operator_Scan_105_Fano.py : Scan of 105 Fano families; confirms only ID-69 satisfies the factorizations
- Bridging_Hseq_HFano_MIS.py : Multiple importance sampling bridging between H_seq and H_Fano
- Fano_Graphon_Verification_Suite.py : Cut-norm convergence, spectral ratio Λ1/Λ0 ≈ 1.99765

## Reproducibility
Fixed random seed: 20260522 (reported in each script)
Numerical precision: 50 decimal digits (mpmath)
All scripts are self-contained and require no external data files

 
