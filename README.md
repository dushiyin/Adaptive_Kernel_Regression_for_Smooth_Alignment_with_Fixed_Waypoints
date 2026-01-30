# Adaptive Kernel Regression for Smooth Alignment with Fixed Waypoints

This repository contains the source code for the arXiv preprint:

> **Adaptive Kernel Regression for Constrained Route Alignment:  
Theory and Iterative Data Sharpening**  
> arXiv:2501.01344

The repository implements the **ANW (adaptive Nadaraya-Watson kernel regression estimator)** and  
**DS-ANW (data–sharpened ANW)** methods proposed in the paper.

---

## Source Code Description

This repository contains the source code for the paper  
**"Energy-based segmentation methods for images with non-Gaussian noise"**.

The proposed framework is designed for image segmentation under non-Gaussian
noise conditions and provides adaptive switching between Gaussian and
non-Gaussian models.

---

## Citation

If you use the code in this repository for your research, please cite our paper:

## Abstract

Route alignment design in surveying and transportation engineering frequently involves fixed waypoint constraints, where a path must precisely traverse specific coordinates. While existing literature primarily relies on geometric optimization or control-theoretic spline frameworks, there is a lack of systematic statistical modeling approaches that balance global smoothness with exact point adherence. This paper proposes an Adaptive Nadaraya-Watson (ANW) kernel regression estimator designed to address the fixed waypoint problem. By incorporating waypoint-specific weight tuning parameters, the ANW estimator decouples global smoothing from local constraint satisfaction, avoiding the "jagged" artifacts common in naive local bandwidth-shrinking strategies. To further enhance estimation accuracy, we develop an iterative data sharpening algorithm that systematically reduces bias while maintaining the stability of the kernel framework. We establish the theoretical foundation for the ANW estimator by deriving its asymptotic bias and variance and proving its convergence properties under the internal constraint model. Numerical case studies in 1D and 2D trajectory planning demonstrate that the method effectively balances root mean square error (RMSE) and curvature smoothness. Finally, we validate the practical utility of the framework through empirical applications to railway and highway route planning. In sum, this work provides a stable, theoretically grounded, and computationally efficient solution for complex, constrained alignment design problems.

