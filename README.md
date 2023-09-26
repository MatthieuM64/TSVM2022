# Two-species Vicsek model

Codes used in the scientific publication:<br>
S. Chatterjee, M. Mangeat, C.-U. Woo, H. Rieger, and J. D. Noh, <a href='https://journals.aps.org/pre/abstract/10.1103/PhysRevE.107.024607'>Flocking of two unfriendly species: The two-species Vicsek model</a>, Phys. Rev. E <b>107</b>, 024607 (2023). Preprint available on <a href='https://arxiv.org/abs/2211.10494'>arXiv</a>.

C++ code on numerical simulation of the two-species Vicsek model.<br>
Exportations: dynamics of flocking; time evolution of the order parameters and the mean-square displacement.<br>
Compile: g++ TSVM.cpp -lgsl -lgslcblas -lm -O3 -s -o TSVM.out.<br>
Run: ./TSVM.out -parameter=value.<br>
List of parameters: v0, eta, rho0, LX, LY, tmax, init, RAN (details as comments in the code).<br>
Generate the movie (rectangle geometry 800x100): python figure_TSVM_dynamics.py -parameter=value (parameters: v0, eta, rho0, LX, LY, DT, tmax, init, ran).<br>
Generate the movie (square geometry 512x512): python figure_TSVM_square_dynamics.py -parameter=value (parameters: v0, eta, rho0, LX, LY, DT, tmax, init, ran).
