# timevarying-causal-KGR

Code for implementing a causal inference analysis on a continuous, time varying treatment using my KGR-SKATER framework and G-methods for spatiotemporal causal inference.

Authors: Jeffrey Wu, Alex Franks, Gareth W Peters (UCSB) 
 
In my second project, I will apply our KGR-SKATER procedure to a causal inference analysis. Instead of just modeling respiratory related mortality like in my first project, in this one, I will specifically be evaluating the direct effects of wildfire specific PM 2.5 on respiratory related mortality in California from 2014 to 2019. In order to extend existing causal analyses in this area, the treatment will be defined as a time varying, continuous variable, which adds complications that preclude traditional causal inference methods from being applicable, these nuances are explained in section III of Hernan and Robins' _What If?_ textbook (insert citation). 
