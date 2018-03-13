# cnVpermtest
Epigenetic and mutational landscape in head and neck squamous cell carcinoma


### Introduction

**This Package allows to do permutation test of Copy Number Alteration (CNA) and Super Enhancers (SE):**

1. **First run function _tcgacnv_ or _tcgacnv_ampdel_:** 
    * Be able to do N times permutation of the CNA segments 
        +For each chromosome within each sample, considered the total range of bases sequenced for CNA detection and randomly moved the discrete CNAs within this range;
        +Use the discrete permuted CNAs to find overlaps with SEs;
        +Took the overlap between the “intersection SEs”/SEs in individual ChIP-seq samples and the permuted CNAs, with the test statistic being the total number of common bases;
        +For specific Amplification and Deletion, we have a seperate function _tcgacnv_ampdel_

2. **Second function _permutation_sum_:**
    * Sum the number of overlap bases for each permutation
  
3. **Third function _permhist_:**
    * Make a histogram of the overlap bases for each permutation;
    * Calculated p-values comparing the distribution under the permutation-based null to the observed values

