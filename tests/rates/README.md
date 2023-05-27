This test simulates the unperturbed motion of 11 non-interating particles that reside on 11 streamlines,
with particles_per_streamline=1, and initial eccentricities logarithmically distributed between
e=1e-6 and e=1e-2. 

Particles having the higher e's are in orbits that cross, which would ordinarily cause epi_int_lite to terminate execution
when detected. However a quirk of the orbit-cross-check is that it emits a nonfatal "divide by zero" warning
when there are less than 3 particles per. Which allows epi_int_lite to evolve these particles to completion.

Point is: the orbit-cross-check function is generating an unimportant divide-by-zero message when 
it tries and fails to check whether these non-interacting particles' orbits do cross. 
Which is fine since this test doesnt care whether these orbits cross.
The check_rates.ipynb notebook then confirms that these particles unperturbed orbits evolve at the expected rates due
to the central planet's J2, namely, that notebook shows that the particles' mean anomaly
advances at the expected rate, dM/dt=Kappa, and that d(wt)/dt=Omega-Kappa

To avoid the "divide by zero" warning, I'll revise epi_int_lite not to execute the orbit-cross-check subroutine 
when particles_per_streamline<3, which I'll encode at a later date.
-jmh 
