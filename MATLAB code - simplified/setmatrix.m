function setmatrix
%c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%c  this subroutine calculates the one electron Fock operator
%c --------------------------------------------------------------------
global   rho ...
    nsym nbas mbasis...
    rp rpor  ...
     c ...
    gl gs dgl nskipe ...
    S df_single
Sll(mbasis,mbasis)=0;
Sss(mbasis,mbasis)=0;
PI(mbasis,mbasis)=0;
V(mbasis,mbasis)=0;

for isym=1:nsym
    for ia=1:nbas(isym)
        ja=ia+nskipe(isym);
        ka=ia+2*nskipe(isym);
        index1=ia+nbas(isym)+2*nskipe(isym);
        for ib=1:nbas(isym)
            jb=ib+nskipe(isym);
            kb=ib+2*nskipe(isym);
            index2=ib+nbas(isym)+2*nskipe(isym);

            Sll(ka,kb)=valint(gl(ja,:).*gl(jb,:).*rp);
            Sss(index1,index2)=valint(gs(ja,:).*gs(jb,:).*rp);

            PI(index1,kb)=valint(gs(ja,:).*dgl(jb,:).*rp);
            
            V(ka,kb)=valint(-gl(ja,:).*gl(jb,:).*rpor.*rho);
            V(index1,index2)=valint(-gs(ja,:).*gs(jb,:).*rpor.*rho);
        end
    end
end

S = Sll + Sss;
PI = PI + PI';

df_single = V + c*PI - 2*c*c*Sss;
