# 9 Others

IMPACT-T has a sister program called IMPACT which also does tracking with space charge forces. The difference between the two is that the independent variable in IMPACT is the longitudinal position z while IMPACT-T uses the time t as the independent variable. The advantage of IMPACT-T over IMPACT is that space charge forces need to be evaluated at constant t. This allows IMPACT-T to more accurately model the space charge forces. This difference between IMPACT-T and IMPACT is most noticeable at low energy. That is, with particles near the gun. When there are external fields (solenoid field, RF field, etc.) then potentially IMPACT-T has to do more work to find the external fields for each individual particle location. This makes it slower than the z-based IMPACT code. It should be noted that at the end of the calculation IMPACT-T stops at some time t so that the particles will have different z.

---

# Acknowledgments

We would like to thank Dr. D. Sagan of the Cornell University for making the first LaTeX draft of the user input/output document after his visit to LBNL in 2006. We would like to thank Dr. C. Mitchell for improving the longitudinal CSR wakefield with an integrated Green function method, Drs. D. Mihalcea and P. Piot for the dielectric wakefield, Dr. D. Zheng for helping replace the NR FFT functions with open source FFTpackage, Ms. C-H. Huang for improving collimators. We would also like to thank all of our collaborators and users for comments, suggestions, and feedbacks to improve this code and for exercising this code on real physical problem studies. This work was supported by the U. S. Department of Energy under Contract no. DE-AC02-05CH11231.

---

# References

[1] J. Qiang, S. Lidia, R. D. Ryne, C. Limborg-Deprey, Phys. Rev. Special Topics - Accel. Beams 9, 044204, (2006).

[2] D. H. Dowell and J. F. Schmerge, Phys. Rev. Special Topics - Accel. Beams 12, 074201, (2009).

[3] K. L.F. Bane, "Short-Range Dipole Wakefields in Accelerating Structures for the NLC," SLAC-PUB-9663, 2003.

[4] P. Craievich, T. Weiland, I. Zagorodnov, "The short-range wakefields in the BTW accelerating structure of the ELETTRA linac," ST/M-04/02.

[5] J. Qiang, C. E. Mitchell, R. D. Ryne, NIM-A 682, 49, (2012).

[6] C. E. Mitchell, J. Qiang, R. D. Ryne, NIM-A 715, 119 (2013).

[7] E. L. Saldin, E. A. Schneidmiller, and M. V. Yurkov, Nucl. Instrum. Methods Phys. Res., Sect. A398, 373 (1997).

[8] M. Borland, Phys. Rev. Special Topics - Accel. Beams 4, 070701 (2001).

[9] G. Stupakov and P. Emma, "CSR Wake for a Short Magnet in Ultrarelativistic Limit," SLAC-PUB-9242, 2002.

[10] W. Press, S. Teukolsky, W. Vetterling, and B. Flannery, *Numerical Recipes in FORTRAN*, Cambridge University Press, New York, 1992.

[11] I. Pogorelov, J. Qiang, R. D. Ryne, M. Venturini, A. Zholents, R. Warnock, "Recent Developments in IMPACT and Application to Future Light Sources," in Proceedings of 9th Intern. Comp. Accel. Physics Conf., p. 182, 2006.

[12] P. L. Morton, Particle Dynamics in Linear Accelerators, Ph. D. Thesis, Midwestern Universities Research Association, The Ohio State University, 1963.

[13] G. A. Loew, R. H. Miller, R. A. Early and K. L. Bane, SLAC-PUB-2295 (1979).

[14] H. A. Enge, Rev. of Sci. Instr. 35, 278, (1964).
