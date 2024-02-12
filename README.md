# Albeck Lab's MATLAB implementation of ARCOS.

This implementation is based on the code presented in https://rupress.org/jcb/article/222/10/e202207048/276138/Automatic-detection-of-spatio-temporal-signaling and available on https://github.com/dmattek/ARCOS with modifications to improve "active" cell detection. 

Our ARCOS-MATLAB port is used to measure cytokine-mediated SPREAD events in human airway epithelial cells in DeCuzzi et al (https://www.biorxiv.org/content/10.1101/2024.02.03.578773v1).

We wish to thank the Pertz Lab and the team that developed ARCOS for this open-source software and allowing us to make a version in MATLAB.

We developed a MATLAB implementation of Automated Recognition of Collective Signaling (ARCOS), which was developed by the Pertz Lab. Following the ARCOS methodology, ERK pulses detected in time series data were binarized and subsequently clustered with Density Based Spatial Clustering of Applications with Noise (DBSCAN). Cluster labels were tracked using the k-nearest neighbors algorithm to match neighbors of clusters within a given distance between time points (epsilon).

Our implementation expands on ARCOS by using pulse detection methods to identify active ERK pulses for binarization, adding automatic estimations of epsilon and the minimum number of points within epsilon distance of a point for that point to be considered a core point, based on the distributions of the input data and cluster area calculations using MATLAB’s native boundary tracing method. Cluster data can be filtered by duration, start time, maximum area, maximum count, mean rate of change of count and mean rate of change of area, effectively narrowing down the pool of spread candidates. We validated our approach by generating synthetic spreads and using Adjusted Rand Indices (ARI) to compare per-timepoint ground truth cluster assignments to ARCOS cluster assignments of varying lifetime and frequency, observing mean per-timepoint ARI values between 0.96 and 0.99.
