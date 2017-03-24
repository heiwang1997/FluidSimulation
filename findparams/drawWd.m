V_PA = 27;
V_PB = 1;
V_RTM = 7.2;
r = 0.05:0.01:0.75;
y = -2 * V_PA + (V_PB * V_PB * V_RTM) ./ (r .* (V_PB - r) .* (V_PB - r));
plot(r, y);

