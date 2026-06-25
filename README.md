# Dynamic-Analysis-of-a-Ground-Vehicle

## Overview
This project investigates the dynamic behavior of a ground vehicle using a combination of finite‑element analysis (HyperMesh-OptiStruct) and numerical modeling in MATLAB. The workflow includes modal, static, and transient analyses in HyperMesh, followed by the construction of reduced‑order models, modal analysis, and numerical time integration using Newmark and Central Difference methods.

## Part 1 — HyperMesh FEA (Summary Only)
- CAD import and refinement
- Modal, static, and transient (modal and direct) analyses
- Export of modal transient response data in selected points (B1-B4 and I1-I3)

## Part 2 — Numerical Model (Included in Repository)
- Assembly of global reduced-order M, K, C matrices
- Integration of 4 suspension-wheel subsystems
- Modal analysis (first 15 natural frequencies)
- Comparison with HyperMesh modal analysis results
- Newmark and Central Difference numerical integration
- Comparison with HyperMesh modal transient results

## Repository Contents
- Modal_Analysis.m
- Transient_Response_Newmark.m
- Transient_Response_Central_Difference.m
- data/ (HyperMesh‑exported matrices and comparison data)

## How to Run
- Load the reduced matrices from the data/ directory
- Run Modal_Analysis.m to compute and compare natural frequencies
- Run Transient_Response_Newmark.m or Transient_Response_Central_Difference.m to simulate the dynamic response
- Compare numerical results with the modal transient FEA data
