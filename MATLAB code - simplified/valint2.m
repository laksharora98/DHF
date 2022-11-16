
function result = valint2(v,f,g)
% This calculates the inner integral in double integrals
global n r rp

result(n) = 0;
u(n) = 0;
for i = 1:n
    r1 = r(i);
    for j =1:n
        r2 = r(j);
        if (j<i)
            u(j) = r2^v/r1^(v+1);
        elseif (j==i)
            u(j) = 1/r1;
        else
            u(j) = r1^v/r2^(v+1);
        end
    end
    result(i) = valint(u.*f.*g.*rp);
end