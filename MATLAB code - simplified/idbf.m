%c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%c function idbf
%c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%c     calculates factorial
%c-----------------------------------------------------------------------
function fac = idbf(x,i)
k = 1;
for m = x
    if(m<-1)
        fprintf(stdout,'factorial problem: m<-1 %d',m);
        return

    else

        if(m==0 || m==-1)

            if(m==0)
                tmp=1;
            end

            if(m==-1)
                if(i==2)
                    tmp=1;
                else
                    fprintf('factorial problem: m<0 %d',m);
                    return
                end
            end

        else

            tmp=1;

            for j=m:-i:1
                tmp=j*tmp;
            end

        end

    end
fac(k) = tmp;
k = k+1;
end
end