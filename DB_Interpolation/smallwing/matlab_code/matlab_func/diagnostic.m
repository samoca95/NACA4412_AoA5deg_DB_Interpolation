function [d99,Uinf,i99,H12,d1,d2] = diagnostic(y,U,u)

if y(1)==0
    y=y(2:end);
    U=U(2:end);
    u=u(2:end);
end
uU=abs(u./U);
lev=0.02;%65;% interp1(U/Uinf0/.99,uU/sqrt(H120),.99)

h=0;
for go=1:5 % just one step needed apparenty!
    h=h+1;
    if h==1
Uinf=max(U);
d1=trapz([0; y],1-[0; U]/Uinf);
d2=trapz([0; y],([0; U]/Uinf).*(1-[0; U]/Uinf));
H12=d1/d2;

%close all
% plot(U/Uinf/.99,uU/sqrt(H12)); hold on
%xlim([0 1])
    else

i99=find(uU/sqrt(H12)<lev,1,'first');

% MA:
if(i99<5)
    buff=find(uU/sqrt(H12)<lev); 
    i99=buff(2); 
end

Uinf=interp1(uU(1:i99+1)/sqrt(H12),U(1:i99+1),lev)/0.99;

% 
% MA:
%disp(max(U)/Uinf)
%disp(d2)


d99=interp1(uU(1:i99+1)/sqrt(H12),y(1:i99+1),lev);
d1=trapz([0; y(1:i99)],1-[0; U(1:i99)]/Uinf);

% disp(min((1-[0; U(1:i99)]/Uinf)))
% disp(min([0; U(1:i99)]))
% disp(min(([0; U(1:i99)]/Uinf)))
% disp(min(([0; U(1:i99)]/Uinf).*(1-[0; U(1:i99)]/Uinf)))

d2=trapz([0; y(1:i99)],([0; U(1:i99)]/Uinf).*(1-[0; U(1:i99)]/Uinf));
H12=d1/d2;
Red2=Uinf*d2;
Red1=Uinf*d1;
%loop
%plot(U/Uinf,uU/sqrt(H12),'r'); 
    end
end