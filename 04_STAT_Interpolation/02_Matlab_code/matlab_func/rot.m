function [xprofu,yprofu,xprofl,yprofl]=rot(xu,yu,xl,yl,yn,alphau,alphal)

for i=1:length(xu)
    i
    Ralpha=[cos(alphau(i)) -sin(alphau(i)); sin(alphau(i)) cos(alphau(i))];
    Ralphapp=[cos(alphal(i)) sin(alphal(i)); sin(alphal(i)) -cos(alphal(i))];
    for j=1:length(yn)
        prod=Ralpha*[0;yn(j)];
        xprofu(j,i)=prod(1,1)+xu(i);
        yprofu(j,i)=prod(2,1)+yu(i);
        prod=Ralphapp*[0;yn(j)];
        xprofl(j,i)=prod(1,1)+xl(i);
        yprofl(j,i)=prod(2,1)+yl(i);
    end
end

end