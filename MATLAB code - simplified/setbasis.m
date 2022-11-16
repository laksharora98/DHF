function setbasis
% This subroutine calculates the GTOs
global stdout nsym nbas nocorb alpha1 beta1 orbj kap r nbasis gl gs dgl ...
    dg1 dg2 kappe orbje kappc orbjc iqe iqc lorba lorbb nskipe nskipc alpha

lorba=[0 1 1 2 2 3 3 4 4 5 5 6 6 7 7];
lorbb=[1 0 2 1 3 2 4 3 5 4 6 5 7 6 8];
nk = lorba+1;
nskipe = [0 cumsum(nbas)];

for isym=1:nsym
    ii = nskipe(isym);
    nn = nbas(isym);
    if (nn~=0)
        kappe(ii+1:ii+nn)=kap(isym);
        orbje(ii+1:ii+nn)=orbj(isym);
        nke(ii+1:ii+nn)=nk(isym);
        iqe(ii+1:ii+nn) = (-1)^(isym+1);
        alpha(ii+1:ii+nn) = alpha1.*(beta1.^(0:(nbas(isym)-1)));
    end
end

nskipc = [0 cumsum(nocorb)];

for isym=1:nsym
    ii = nskipc(isym);
    nn=nocorb(isym);
    if(nn~=0)
        kappc(ii+1:ii+nn)=kap(isym);
        orbjc(ii+1:ii+nn)=orbj(isym);
        iqc(ii+1:ii+nn) = (-1)^(isym+1);
    end
end

fprintf(stdout,' total number of basis function= %d',nbasis);
fprintf(stdout,'symmetry   no. of occ orbits   nskipc');

for isym=1:nsym
    fprintf(stdout,'%d            %d            %d',isym,nocorb(isym),nskipc(isym));
end

fprintf(stdout,'symmetry   no. of basis function   nskipe');
for isym=1:nsym
    fprintf(stdout,'%d            %d            %d',isym,nbas(isym),nskipe(isym));
end

ifact=idbf(2.*nke-1,2);
power=(2.*nke+1)./2;
cnl=sqrt((2.^(2.*nke+1.5)).*(alpha.^power)./(ifact.*sqrt(pi)));
fact1=4.*alpha.*((nke+kappe).^2)./(2.*nke-1);
fact2=(2.*nke+1).*alpha;
fact3=-4.*alpha.*(nke+kappe);
cns=sqrt(1./(fact1+fact2+fact3));

expon= exp(-(alpha.').*r.*r);

gl=(cnl.').*expon.*(r.^(nke.'));
gs=(cns.').*gl.*((nke.'+kappe.')./r - 2.*(alpha.').*r);
gs(:,1)=0;
dgl=gl.*((nke.'+kappe.')./r - 2.*(alpha.').*r);
dgl(:,1)=0;
dg1=((nke.')./r) - (2.*(alpha.').*r);
dg1(:,1)=0;
dg2=-(((nke.'+kappe.')./(r.^2))+2.*(alpha.'))./...
    ((nke.'+kappe.')./r-2.*(alpha.').*r);
dg2(:,1)=0;

end