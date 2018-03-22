# Self-avoiding-walk-of-DNA
Two complementary single strands of DNA are considered as self avoiding random walkers. And each nitrogenous base of the DNA strand is considered as the link of the random walker. Each strand is a random walker with a step size of DNA length (length_dna), and the bases are represented as "sdna_one" and "sdna_two". Whenever the nitrogenous bases of the complementary strands occupy same position it is considered as a "contact". If the contacted bases are non complementary then the contact is called as non correct contact (no_corr_cont) and if the contacted bases complementary then it is called correct contact (corr_cont).

## Variations
The fraction of correct contacts can be calculated for a range of DNA length (bps_strart and bps_end), for example in this code I have calculated for DNA strands of length 5 to 200. This helps to correlate the number of correct contacts or the nucleus size at zeroth time with the DNA length. 
Similarly, the confinement volume as well as the dimension (2 or 3) can be changed. The measures of the confinement can be given by min and max. Through this the effect of confinement volume and the dimension on the nucleus size can be calculated. 
By simulating the DNA strands for a particular measurement for 10^5 times, this code will cover significant number of random walker conformations.


## Authors

* **Niranjani** 


## Acknowledgments

* Dr.R.Murugan, IITM - The core concept of the random walk was developed by him.
