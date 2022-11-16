%c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%c subroutine scfiter
%c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%c     this performs the scf iteration.
%c------------------------------------------------------------------
function scfiter

global stdout     ...
    tiny      ...
    nsym nbas nocorb  delta  ...
    crit maxit npower orbj  ...
    mbasis ...
    nskipe  ...
    S df_single df_two ...
    en sa drs vec

nh_i = [' S ', ' P-', ' P ', ' D-', ' D ', ' F-', ' F ' ,' G-' ,' G ', ' H-', ' H ', ...
    ' I-', ' I ', ' K-',  ' K '];
oe(mbasis,nsym)=0;

iter=0;
p1(mbasis,mbasis)=0;

fprintf(stdout,['convergence data\n maximum no. of iteration= %6d ...' ...
    'convergence criterion   = 1.0d-%2d\n'],maxit,npower);
fprintf(stdout,'cycle      density conv\n');

while(iter<maxit)
    iter=iter+1;
    if(iter~=1)
        dfpotiter();
    end
    for isym=1:nsym
        ndim=nbas(isym);
        ndim2=2*ndim;
        skip = 2*nskipe(isym);
        stemp = S(1+skip:ndim2+skip,1+skip:ndim2+skip);

        drs = stemp;

        [drs,temp] = eig(drs);
        [sa,ind] = sort(diag(temp),'descend');
        drs = drs(:,ind);

        %c-------------------------------------------------------------------
        %c       canonical transformation
        %c------------------------------------------------------------------
        
        if (min(abs(sa))>tiny)
            can = drs./sqrt(abs(sa'))
        else
            fprintf('there is diagonalisation problem %f\n\n',min(abs(sa)));
            return
        end

        if(iter~=1)
            a=df_single(1+skip:ndim2+skip,1+skip:ndim2+skip)+...
                df_two(1+skip:ndim2+skip,1+skip:ndim2+skip)
        else
            a=df_single(1+skip:ndim2+skip,1+skip:ndim2+skip)
        end

        df_mat = can'*a*can

        [df_mat,temp] = eig(df_mat)
        [eng,ind] = sort(diag(temp),'ascend')
        df_mat = df_mat(:,ind)

        drs = df_mat

        eigv = can*drs

        oe(1:ndim2,isym)=eng 

        pold(1+skip:ndim2+skip,1+skip:ndim2+skip) = p1(1+skip:ndim2+skip,1+skip:ndim2+skip);

        nn=nocorb(isym);
        if(nn~=0)
            p1(1+skip:ndim2+skip,1+skip:ndim2+skip) = (2*orbj(isym)+1)*eigv(:,1:nn)*eigv(:,1:nn)'
        end

        vec(1+skip:ndim2+skip,1+skip:ndim2+skip)=eigv

    end
    delta=sum(sum((p1-pold).^2))

    delta=sqrt(delta/4.0)

    fprintf(stdout,'%3d   %20.7f\n',iter,delta);
    fprintf('%3d   %20.7f\n',iter,delta);
    if(delta<crit)
        break
    end
end
if (iter>=maxit)
    fprintf(stdout,'\n\n scf fails to converges at cycle %4d\n',iter);
    return
end

en = 0.5*sum(sum(p1.*(2*df_single+df_two)))

fprintf(stdout,'\n\nscf converges at cycle, %4d\n',iter);
fprintf(stdout,'\n\nscf electronic energy=, %25.15f\n',en);

fprintf(stdout,'\n\norbital energies (+ve and -ve)\n');

for isym=1:nsym
    for ibas=1:nbas(jsym)
        fprintf(stdout,['%2d %s %25.15f %25.15f\n',ibas,nh_i(isym), ...
            oe(ibas,isym),oe(ibas+nbas(jsym),isym)]);
    end
end

end