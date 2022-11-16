function setgrd
%cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%cc  this subroutine defines the grids                           cc
%cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

global stdout rp rpor rnt h n r

eph = exp(h);

%c ---------------------------------------------------------------------
%c   set up the arrays r, rp, rpor
%c ---------------------------------------------------------------------

ett = eph.^(0:n-1);
ettm1 = ett - 1;
r = rnt*ettm1;
rp = rnt*ett;
rpor = ett./ettm1;
rpor(1) = 0;

fprintf(stdout,'stepsize=%f total grids= %d\n',h,n);
fprintf(stdout,'maximum radius r(n) = %f\n',r(n));

end

