function [res] = trapezoidal2(x,y,v,ang,ind)

for i = 1:length(x)-1
    dx(i) = x(i+1)-x(i);
    dy(i) = y(i+1)-y(i);
    dv(i) = (v(i+1)+v(i))/2;
    dang(i) = (ang(i+1)+ang(i))/2;
end

if ind == 1
    res = sum(dx.*dv);
elseif ind == 2
    res = sum(dx.*dv.*tan(dang));
    % res = sum(dy.*dv);
elseif ind == 3
    % res = sum(dy.*dv./tan(dang));
    res = sum(dy.*dv);
elseif ind == 4
    res = sum(dy.*dv./tan(dang));
    % res = sum(dx.*dv);
elseif ind == 5 
    res = sum(dx.*dv.*(dx/2+x(1:end-1)));
elseif ind == 6 
    res = sum(dx.*dv.*(dy/2+y(1:end-1)));
elseif ind == 7 
    res = sum(dx.*dv.*tan(dang).*(dx/2+x(1:end-1)));
else 
    res = sum(dx.*dv.*tan(dang).*(dy/2+y(1:end-1)));
end