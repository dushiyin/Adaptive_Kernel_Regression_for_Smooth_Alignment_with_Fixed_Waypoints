# Adaptive_Kernel_Regression_for_Smooth_Alignment_with_Fixed_Waypoints

## Source Code for "Adaptive Kernel Regression for Constrained Route Alignment"

This repository contains the source code for the arXiv preprint:

> **Adaptive Kernel Regression for Constrained Route Alignment: Theory and Iterative Data Sharpening**  
> arXiv:2501.01344

The repository implements the **ANW (adaptive Nadaraya-Watson kernel regression estimator)** and **DS-ANW (data–sharpened ANW)** methods proposed in the paper.


## Citation

If you use the code in this repository for your research, please cite our paper:

S. Du, Y. Chen, W. Yang, Q. Li, and X. Shi, Adaptive kernel regression for constrained route alignment: Theory and iterative data sharpening, arXiv preprint arXiv:2601.01344, (2026).

## Abstract

Route alignment design in surveying and transportation engineering frequently involves fixed waypoint constraints, where a path must precisely traverse specific coordinates. While existing literature primarily relies on geometric optimization or control-theoretic spline frameworks, there is a lack of systematic statistical modeling approaches that balance global smoothness with exact point adherence. This paper proposes an Adaptive Nadaraya-Watson (ANW) kernel regression estimator designed to address the fixed waypoint problem. By incorporating waypoint-specific weight tuning parameters, the ANW estimator decouples global smoothing from local constraint satisfaction, avoiding the "jagged" artifacts common in naive local bandwidth-shrinking strategies. To further enhance estimation accuracy, we develop an iterative data sharpening algorithm that systematically reduces bias while maintaining the stability of the kernel framework. We establish the theoretical foundation for the ANW estimator by deriving its asymptotic bias and variance and proving its convergence properties under the internal constraint model. Numerical case studies in 1D and 2D trajectory planning demonstrate that the method effectively balances root mean square error (RMSE) and curvature smoothness. Finally, we validate the practical utility of the framework through empirical applications to railway and highway route planning. In sum, this work provides a stable, theoretically grounded, and computationally efficient solution for complex, constrained alignment design problems.

## File Structure

This repository contains the complete implementation of the proposed Adaptive Nadaraya–Watson (ANW) regression and iterative data-sharpening methods, along with baseline methods used for comparison in both simulation studies and real-world applications.


### `application/`

This folder contains scripts for **real-world route alignment applications**.

- `high-speed_rail.R`  
  Implements the proposed methods on high-speed rail alignment data.

- `highway_all_methods.R`  
  Applies all competing methods (including NW and ANW-based approaches)
  to highway alignment data for comprehensive comparison.

- `highway_all_methods_map.R`  
  Generates map-based visualizations of the highway alignment results.

- `highway_anw_vs_nw.R`  
  Focuses on the comparison between the standard Nadaraya–Watson estimator
  and the proposed Adaptive Nadaraya–Watson method for highway data.

- `rail_all_methods.R`  
  Applies all methods to railway alignment data.

- `rail_anw_vs_nw.R`  
  Direct comparison between NW and ANW methods on railway alignment data.

---

### `simulation/`

This folder contains scripts for **controlled simulation studies** designed
to evaluate the numerical properties of the proposed methods.

- `1D_all_methods.R`  
  One-dimensional simulation comparing all competing methods.

- `1D_anw_vs_nw.R`  
  One-dimensional simulation focusing on NW versus ANW comparison.

- `2D_all_methods.R`  
  Two-dimensional simulation comparing all competing methods.

- `2D_anw_vs_nw.R`  
  Two-dimensional simulation focusing on NW versus ANW comparison.


## Usage Instructions

## Contact

If you have any questions, please contact us at [dusy77@student.ubc.ca].

