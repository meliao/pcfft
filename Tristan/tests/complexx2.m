
function [r,d,d2] = complexx2(t,cpars)
%complexx
% return position, first and second derivatives of parameterized
% complexified layered media interface on the x-axis from [-L,b]
% or [b,L].
%
% Inputs:
% t - parameter values to evaluate these quantities
%
% Optional inputs:
% r - coordinates
% d - first derivatives w.r.t. t
% d2 - second derivatives w.r.t. t

L = cpars.L; c1=cpars.c1; c2=cpars.c2;

% if icurve==1 % left part
%     xs = t-c1*1i*erfc((t+L)/c2);
%     xp = 1+c1*1i*2/sqrt(pi)/c2*exp(-(t+L).^2/c2^2);
%     xpp = -c1*1i*4/sqrt(pi)/c2^3*(t+L).*exp(-(t+L).^2/c2^2);
% elseif icurve==2 % right part
%     xs = t+c1*1i*erfc((L-t)/c2);
%     xp = 1+c1*1i*2/sqrt(pi)/c2*exp(-(L-t).^2/c2^2);
%     xpp = -c1*1i*4/sqrt(pi)/c2^3*(t-L).*exp(-(L-t).^2/c2^2);   
% end

% erfci = @(t) t.*erfc(t) - exp(-t.^2)/sqrt(pi) - 2*t;
% xs = t-c1*1i*erfci((t+L)/c2)+c1*1i*erfci((L-t)/c2);
% xp = 1-c1*1i*erfc((t+L)/c2)+c1*1i*erfc((L-t)/c2);
% xpp = c1*1i*2/sqrt(pi)/c2*exp(-(t+L).^2/c2^2)+c1*1i*2/sqrt(pi)/c2*exp(-(L-t).^2/c2^2);
% xpp = -c1*1i*4/sqrt(pi)/c2^3*(t+L).*exp(-(t+L).^2/c2^2)-c1*1i*4/sqrt(pi)/c2^3*(t-L).*exp(-(L-t).^2/c2^2);   



xs = t-c1*1i*erfc((t+L)/c2)+c1*1i*erfc((L-t)/c2);
xp = 1+c1*1i*2/sqrt(pi)/c2*exp(-(t+L).^2/c2^2)+c1*1i*2/sqrt(pi)/c2*exp(-(L-t).^2/c2^2);
xpp = -c1*1i*4/sqrt(pi)/c2^3*(t+L).*exp(-(t+L).^2/c2^2)-c1*1i*4/sqrt(pi)/c2^3*(t-L).*exp(-(L-t).^2/c2^2);   


ys = zeros(size(t));
yp = ys;
ypp = ys;

r = [(xs(:)).'; (ys(:)).'];
d = [(xp(:)).'; (yp(:)).'];
d2 = [(xpp(:)).'; (ypp(:)).'];

end

