# MIT_BIH_QRS_ECG_DETECTION

This package includes all the code to generate the results in the following paper:

D. McCubbins, B. R. Thapa, K. Hall, and J. Bae. "Simple and Effective Signal Processing Techniques to Detect Heartbeats from Electrocardiogram." 2025 IEEE SoutheastCon (in press).

If you use any of the code published in this repository, partially or fully, you MUST cite the above conference paper.

The primary script in this repository is:

- **`ECGDetection.m`**  
  This script implements the complete signal processing pipeline for detecting R-peaks from raw ECG signals. The pipeline includes: 
  - **Band-pass filtering**  
  - **Peak detection and QRS complex identification**  
  - **QRS duration measurement**  
  - **Wide Complex Tachycardia (WCT) / Narrow Complex Tachycardia (NCT) determination**
 
  The code is implemented in MATLAB and is designed for clarity, reproducibility, and ease of use with standard MATLAB toolboxes.


Data files (ECG) can be found at the following website:
https://www.kaggle.com/datasets/protobioengineering/mit-bih-arrhythmia-database-modern-2023 

The original MIT-BIH Arrhythmia Database can be found at the following website: 
https://www.physionet.org/content/mitdb/1.0.0/

A proper citation is required to use the data set, as described in the website.
