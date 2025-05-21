function [dsfx,sfx]=splinediff6(fx,x,step)

sfx=spline(x(1:step:end),fx(1:step:end),x);
dsfx=diff6(sfx,x);

end