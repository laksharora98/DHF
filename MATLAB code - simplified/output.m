function output()

global stdout iwfn ...
    mns rho ...
    nsym nbas nocorb ...
r ...
    rp rpor h n ...
    nskipe ...
    eng1 en

icntb = zeros(mns);
icntc = zeros(mns);

fprintf(stdout,['  these are the orthogonality of the orbitalsn\n\nsl. no. ...' ...
    '   <a|b>          <a|r|b>         <a|r*r|b>        <a|(1/r)|b>\n\n']);

fprintf(iwfn,'%f %d\n',h,n);
fprintf(iwfn,'%f %f %f\n',r,rp,rpor);

%c     write only the relevant orbitals i.e skip writting the
%c     negavtive energy orbitals

for isym=1:nsym

    icntb(isym)=0;
    icntc(isym)=0;

    for ia=1:nbas(isym)
        ka=ia+nskipe(isym);

        if(temp2)
            icntb(isym)=icntb(isym)+1;
            fprintf(fort15,'%d %f\n',ia,eng1(ka));
            fprintf(iwfn,'%f\n',eng1(ka));
            fprintf(iwfn,'%d %d\n',pf(1:n,ka),qf(1:n,ka));
        end
    end
end

fprintf('%20s %4d','no. of sym:\n', nsym);
    fprintf('%20s %4d','no. of basis:\n',icntb(1:nsym));
fprintf('%20s %4d','no. of occ. orbs:\n',nocorb(1:nsym));
fprintf(fort15,'%20s %4d','no. of sym:\n', nsym);
    fprintf(fort15,'%20s %4d', 'no. of basis:\n',icntb(1:nsym));
fprintf(fort15,'%20s %4d','no. of occ. orbs:\n',nocorb(1:nsym));
    fprintf(iwfn,'%15.6f\n',rho);
end