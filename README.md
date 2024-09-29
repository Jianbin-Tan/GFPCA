# Graphical Functional Principal Component Analysis (GFPCA)

This README accompanies the paper "Graphical Principal Component Analysis of Multivariate Functional Time Series" authored by Jianbin Tan, Decai Liang, Yongtao Guan, and Hui Huang. The paper is accessible at [Taylor & Francis Online](https://www.tandfonline.com/doi/full/10.1080/01621459.2024.2302198).

## Citation
If you use this data or methodology in your work, please cite:

```bibtex
@article{tan2024graphical,
  title={Graphical principal component analysis of multivariate functional time series},
  author={Tan, Jianbin and Liang, Decai and Guan, Yongtao and Huang, Hui},
  journal={Journal of the American Statistical Association},
  pages={1--24},
  year={2024},
  publisher={Taylor & Francis}
}

## 1. Data
### Abstract

This dataset comprises hourly PM2.5 concentration readings collected from 24 monitoring stations across three cities in China during the winter of 2016. To stabilize variance, a square-root transformation was applied to the data. The study spans a total of 60 days.

### Availability
All data required to replicate our findings are available.

### Data Dictionary
The dataset, stored in "fda_dat.rda" within the "Data" directory, includes:
- Hourly PM2.5 concentration measurements from 24 monitoring stations over 60 days.
- Descriptive details about each monitoring station, including station name, coordinates (longitude and latitude), province, and city.

## 2. Code
### Abstract
The primary codebase for this paper is provided in "Function.R". Additionally, a detailed guide, "GFPCA.pdf", explains the functions within "Function.R". An illustrative example demonstrating the use of these functions is available in "Example code.R".

### Reproducibility
- The "GFPCA" function in "Function.R" implements the GFPCA method described in this paper.
- Curve reconstruction procedures are detailed in "Example code.R".

## 3. Proof
### Abstract
The proofs supporting the methodologies used in this paper are detailed in "Supplementary Materials of Graphical Principal Component Analysis of Multivariate Functional Time Series.pdf".