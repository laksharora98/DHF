% Input Values
%
%   (a) atom : symbol of atom in two capital letters
%
%   (b) nsym:  no. of symmetry to be used
%
%   (c) nbas(i) :  no. of basis for each symmetry
%
%   (d) nocorb(i) : no. of occupied orbitals for each symmetry
%
%   (e) alpha and beta : values for generating basis
%
%   (f) maxit : max. no. of iterations it should check for convergence
%
%   (g) npower : tolerence parameter for 10^(-npower)
%
%   (h) amass : mass number
%
%   (i) z : atomic number

%---------------------------------------------------------------------
% the output files are:                                               
%                                                                     
%     scfout.out:  contains all information of orbital energies       
%     wfn.dat  :   binary file containing df wave functions           
%                                                                     
%---------------------------------------------------------------------


% this is the tiny number to check density consistency in each iteration
global tiny
tiny = 1.0e-12;

readinp();
setgrd();
setbasis();
setmatrix();
scfiter();

