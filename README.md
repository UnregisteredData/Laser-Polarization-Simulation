# Laser Polarization Simulation

This code takes inputs of Chi-e, and initial polarization, although Phi and Phi-CE values can also be inputted, if desired.

This code works by splitting up the integral found in eq. 1 into four different integrals. Because computers solve integrals numerically, there were issues withthe integral being divergent or slowly convergent for certain values of Chi-e. This was able to be avoided with a series of substitutions. 

Although values of Phi, Phi-CE (pulse duration and carrier envelope phase, respectively) can be manually inputted, previous research by (Seipt et al., 2018) have found that a the asymmetry of the magnetic field (~0.17) is maximal when Phi-CE = pi/2, and Phi = 3.27. Thus, in our simulations, we did not alter these values. 

Feel free to test the code out and play with different values of Chi-e, initial polarization, Phi, and Phi-CE.

Be warned, it still gives a divergent/slowly convergent warning, but still spitsout a value. It still takes a few seconds to generate values though. :)
