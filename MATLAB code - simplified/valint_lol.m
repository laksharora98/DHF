function res = valint_lol(fx)
%c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%c    this is the integration routine using numerical formula       *
%c simpson's formular truncate at h^4, but here it is at h^13       *
%c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
global h n


res=0;

for i=n+1:n+1
    fx(i)=0;
end

for i=1:2:n
    res=res+(fx(i)+fx(i+1));
end

res=h*(res/2);

end
