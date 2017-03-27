V_PA = 27;
V_PB = 1;
V_RTM = 7.2;
r = 0.05:0.01:0.75;
%r = [0.11685700197446966 0.5410226704894886]';
% y = -2 * V_PA + (V_PB * V_PB * V_RTM) ./ (r .* (V_PB - r) .* (V_PB - r));
y = -2 * V_PA .* r + V_RTM * log(r ./ (V_PB - r)) + V_RTM * V_PB ./ (V_PB - r);
figure;
plot(r, y);

