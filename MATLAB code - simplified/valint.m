function res = valint(fx)
%c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%c    this is the integration routine using 11-point Newton-Cotes
%c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
global h n


res=0;

for i=n+1:n+10
    fx(i)=0;
end

for i=1:10:n
    res=res+(16067*(fx(i)+fx(i+10))+427368*fx(i+5) ...
        +106300*(fx(i+1)+fx(i+9))-48525*(fx(i+2)+fx(i+8)) ...
        +272400*(fx(i+3)+fx(i+7))-260550*(fx(i+4)+fx(i+6)));
end

res=5.0*h*(res/299376.0);

end
