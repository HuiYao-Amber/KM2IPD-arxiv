# Benchmarking IPD reconstruction from Kaplan–Meier curves

This repository contains the R simulations for comparing four methods that reconstruct individual participant data (IPD) from Kaplan–Meier curves: iKM, IPDfromKM, KMtoIPD, and QP. The single-trial benchmark measures reconstruction error. The meta-analysis benchmark evaluates how substituting reconstructed estimates for unavailable reported hazard ratios affects pooled estimates.

## Run a small example

The scripts were tested with R 4.5.2. Run these commands from the repository root. The installer checks the versions of the direct dependencies and installs missing packages in a local `.Rlib/` directory.

```sh
Rscript scripts/install_dependencies.R
Rscript tests/smoke_test.R
Rscript scripts/run_simulations.R --part=single
Rscript scripts/run_simulations.R --part=meta
```

The single-trial example runs all four methods on the same simulated trial. The meta-analysis example simulates three trials and evaluates two proportions of unavailable hazard ratios. Results are written under `outputs/`. Use `--output=PATH` to choose a new directory for another run. The smoke test checks all four methods and both simulation workflows.

## Run the main benchmarks

The R code generates the study designs and simulated inputs. [Benchmark settings](configs/benchmark.R) define the method options and seed bases. The following commands use the replicate counts from the study.

```sh
# 61 single-trial settings, 1,000 replicates per setting
Rscript scripts/run_simulations.R --part=single --full=true --reps=1000

# 45 main and 16 stress scenarios, 300 replicates per scenario
Rscript scripts/run_simulations.R --part=meta --full=true --preset=main --reps=300 --output=outputs/meta_main
Rscript scripts/run_simulations.R --part=meta --full=true --preset=stress --reps=300 --output=outputs/meta_stress
```

The single-trial results record reconstruction and evaluation success separately. Keep failed attempts when calculating success rates. Meta-analysis outputs include trial-level estimates and pooled estimates by strategy. The full simulations take substantially longer than the examples.

## Sensitivity analysis

```sh
Rscript scripts/run_sensitivity.R --part=single-factor
Rscript scripts/run_sensitivity.R --part=input-ablation
```

These commands run illustrative designs. The first varies digitization error and the second varies the inputs available to reconstruction methods. Both write the design and the results to the output directory. The second also writes the variants and their settings. These examples show the input columns and output format. Reconstruction failures are recorded in the results.

For a full sensitivity run, provide your own design with `--full=true --design=PATH`. The input-ablation workflow also requires `--qp-variants=PATH` specifying the QP variants. Reproducing the paper's complete sensitivity analyses requires the corresponding study inputs, which are not distributed in this release.

## Method implementations

| Workflow | Benchmark implementation | Inputs used in the main comparison |
| --- | --- | --- |
| iKM | [reconstructKM adapter](src/methods/ikm.R) | KM coordinates and numbers at risk |
| IPDfromKM | [Single-trial reconstruction](src/simulation/single_trial.R) | KM coordinates, numbers at risk, and total events by arm |
| KMtoIPD | [Single-trial reconstruction](src/simulation/single_trial.R) | KM drop coordinates, initial arm size, and censor marks |
| QP | [CIFresolve adapter](src/methods/qp.R) | KM coordinates, numbers at risk, and censor marks |

The [single-trial module](src/simulation/single_trial.R) simulates event and censoring times and constructs digitized inputs. The [meta-analysis module](src/simulation/meta_analysis.R) generates trial collections and pools estimates. The [sensitivity module](src/simulation/sensitivity.R) runs comparisons with different input settings. All participant records are simulated locally.
