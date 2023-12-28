#  Graphical Functional Principal Component Analysis (GFPCA)

This file contains the data, codes, and proof from the paper "Graphical Principal Component Analysis of Multivariate Functional Time Series" by Jianbin Tan, Decai Liang, Yongtao Guan, and Hui Huang.


------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## 1.  Data
### 1) Abstract

We collected hourly PM2.5 concentration readings from 24 monitoring stations located in three cities in China during the winter of 2016. The dataset covers a total of 60 days, and we have applied a square-root transformation to stabilize the variance.

### 2) Availability
The data necessary to reproduce our results is available.

### 3) Data dictionary
The dataset is included in "fda_dat.rda" within the "Data" file. It contains hourly measurements of PM2.5 concentration from 24 monitoring stations, covering a period of 60 days. Additionally, the dataset includes information about the monitoring stations, such as station names, longitudes, latitudes, provinces, and cities.

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## 2. Code
### 1) Abstract
We provide the core codes of this article in "Function.R". Besides, we also include a help document named "GFPCA.pdf" for the functions in "Function.R" and provide an example code in "Example code.R" to use the functions.

### 2) Reproducibility
- The R function "GFPCA" in "Function.R" can be used to implement GFPCA in this article. 
- The codes to reconstruct curves are presented in "Example code.R".

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## 3. Proof
### 1) Abstract
We provide the proof of this article in "Supplementary Materials of Graphical Principal Component Analysis of Multivariate Functional Time Series.pdf". 





