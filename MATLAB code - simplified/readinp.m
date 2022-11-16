function readinp
% This subroutine reads input from scfinp.dat
global atom stdout stdin iwfn nsym nbas nocorb alpha0 beta maxit npower ...
    orbj kap amass z nbasis nocc mbasis crit h n rnt rrms rho 

stdout = fopen('scfout.out','w');
stdin = fopen('scfinp.dat','r');
iwfn = fopen('wfn.dat','w');

% read atom
atom = fgetl(stdin); 

% read number of symmetries
nsym = str2double(fgetl(stdin)); 

% read number of basis for each symmetry in the form of space separated values
nbas = str2num(fgetl (stdin)); 

% read number of occupied orbitals for each symmetry in the form of space separated values
nocorb = str2num(fgetl (stdin)); 

% total number of basis
nbasis=sum(nbas); 

mbasis=2*nbasis;

% read constants to generate exponents
alpha0 = str2double(fgetl (stdin)); 
beta = str2double(fgetl (stdin));

% read maximum number of iterations
maxit = str2double(fgetl (stdin));
% read tolerance
npower = str2double(fgetl (stdin));
% read atomic mass
amass = str2double(fgetl (stdin));
% read atomic number
z = str2double(fgetl (stdin)); 
% total number of occupied orbitals
nocc=sum(nocorb); 
% values of j for each symmetry
orbj = [0.5 0.5 1.5 1.5 2.5 2.5 3.5 3.5 4.5 4.5 5.5 5.5 6.5 6.5 7.5]; 
% values of kappa for each symmetry
kap = [-1 1 -2 2 -3 3 -4 4 -5 5 -6 6 -7 7 -8]; 
% number of values in the grid along r
n = 740; 
% value of h to generate exponential grid
h = 3e-2; 
%%h = 739*0.03/(n-1)
% critical value to end iteration procedure
crit=(0.1)^npower; 

cp=2.2677*1e-5*amass^(1/3);
rrms=sqrt(3/5)*cp;

rnt= exp(-65/16)/z;
rho(1:n) = z; % value for point nucleus

fprintf(stdout,'                     ^^^^*                         \n');
fprintf(stdout,'<<<<<<<<<<<<<<<<<<<<<<<<<o>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n');
fprintf(stdout,'                     ^^^^*                         \n');
fprintf(stdout,'scf run for df wave function\n');
fprintf(stdout,atom);
fprintf(stdout,'\n');
fprintf(stdout,'^^^^^^^^^^^^^^^^^^^\n');
fprintf(stdout,'*                                    *\n');
fprintf(stdout,'*                                    *\n');
fprintf(stdout,'*                                    *\n');
fprintf(stdout,'*        relativistic scf run        *\n');
fprintf(stdout,'*                                    *\n');
fprintf(stdout,'*            written by              *\n');
fprintf(stdout,'*                                    *\n');
fprintf(stdout,'*           Laksh Arora              *\n');
fprintf(stdout,'*                                    *\n');
fprintf(stdout,'*                                    *\n');
fprintf(stdout,'*                                    *\n');
fprintf(stdout,'^^^^^^^^^^^^^^^^^^^\n');
fprintf(stdout,' gaussian type of orbitals chosen\n');
fprintf(stdout,' point nucleus is considered          \n');
fprintf(stdout,' neutral atom has been considered     \n');
fprintf(stdout,'\n');
fprintf(stdout,'>>>>>>>>>>>>>>> input data >>>>>>>>>>>>>>>>>>>>>>>>\n');
fprintf(stdout,'number of symmetries = %2d\n',nsym);
fprintf(stdout,'nbas -->');
fprintf(stdout,char(strjoin(string(nbas))));
fprintf(stdout,'\nocc orb  --> ');
fprintf(stdout,char(strjoin(string(nocorb))));
fprintf(stdout,'\n');
fprintf(stdout,'basis type :Universal Basis\n');
fprintf(stdout,'alpha0 = %10.4f beta = %10.4f\n',alpha0,beta);
fprintf(stdout,'max iter = %d conv. = %10.4f\n',maxit, crit);
fprintf(stdout,'atomic number %10.4f and mass = %10.4f\n',z,amass);
fprintf(stdout,'<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n');
fprintf(stdout,'nuclear r0 = %f  rrms = %f',rnt,rrms);
 
end