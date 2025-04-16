function  [x_pts,y_pts,npoints,yn,alphau,alphal,su,sl,x] = create_profiles_wing_interval(m,p,xx,nle,tol,dy1,dy2,yn,r,xstart,xend)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

disp('compute wing profiles')

% Generate a full-resolution NACA profile
% [xubuff, yubuff, alphaubuff, xlbuff, ylbuff, alphalbuff, xbuff] = naca_prof(m, p, xx, nle, tol);
[xubuff, yubuff, alphaubuff, xlbuff, ylbuff, alphalbuff, xbuff] = naca_prof(m, p, xx, round(pi/tol), tol);

% Segment the profile between xstart and xend
buffmask = ((xbuff > xstart) & (xbuff < xend));

x      = xbuff(buffmask);

xu     = xubuff(buffmask);
yu     = yubuff(buffmask);
alphau = alphaubuff(buffmask);
su     = [0, cumsum(sqrt(diff(xu).^2 + diff(yu).^2))];

xl     = xlbuff(buffmask);
yl     = ylbuff(buffmask);
alphal = alphalbuff(buffmask);
sl     = [0, cumsum(sqrt(diff(xl).^2 + diff(yl).^2))];

%Downsample to nle points
x_target = linspace(xstart,xend,nle);
for ix = 1:length(x_target)
    [tmp, closest_index_] = min(abs(xu - x_target(ix)));
    closestValue_ = xu(closest_index_);
    closest_index(ix) = closest_index_;
    closestValue(ix) = closestValue_;
end
x      = x(closest_index);
xu     = xu(closest_index);
yu     = yu(closest_index);
alphau = alphau(closest_index);
su     = su(closest_index);
xl     = xl(closest_index);
yl     = yl(closest_index);
alphal = alphal(closest_index);
sl     = sl(closest_index);

% Store wall-normal profile
ymax = yn;
y1 = 0; % First point is at the wall
% n = ceil(log(-(y1 * (r - 1) + yn * (1 - r) - dy1) / dy1) / log(r) + 1); % Number of points if we were to use r up to yn
% yn(1) = y1;
% for i = 2:n
%     yn(i) = yn(1) + dy1 * (1 - r^(i - 1)) / (1 - r);
% end
yn(1) = y1;
i = 2;
while yn(end) < ymax
    exp_yn = yn(1) + dy1 * (1 - r^(i - 1)) / (1 - r);
    if (exp_yn - yn(i-1)) < dy2
        yn(i) = exp_yn;
    else
        yn(i) = yn(i-1) + dy2;
    end
    i = i + 1;
end

% Rotate and rearrange the interpolating mesh
[xprofu, yprofu, xprofl, yprofl] = rot(xu,yu,xl,yl,yn,alphau,alphal);

npoints = 0;

% Upper surface
for i = 1:length(xprofu(1, :))
    for j = 1:length(xprofu(:, 1))
        npoints = npoints + 1;
        x_pts(npoints, 1) = xprofu(j, i);
        y_pts(npoints, 1) = yprofu(j, i);
    end
end

% Lower surface
for i = 2:length(xprofl(1, :))
    for j = 1:length(xprofl(:, 1))
        npoints = npoints + 1;
        x_pts(npoints, 1) = xprofl(j, i);
        y_pts(npoints, 1) = yprofl(j, i);
    end
end

% Store wall-normal profile yn
for i = 1:length(yn)
    npoints = npoints + 1;
    x_pts(npoints, 1) = -1;
    y_pts(npoints, 1) = yn(i);
end

end

