%% General Physics Model
% Water van der waals constants - JL Notation
a = 1708.95;
b = 590.357;
R = 8.314 / (1.8e-2);

% Normailized form
a = 27;
b = 1;
R = 1;

%% Thermal Control Bar
theta_critical = 8 * a * b / (27 * R);
theta = 0.9 * theta_critical;

%% Fit Curve
rho = 0:0.0001:(0.7 * b);
p = ((R * b * rho * theta) ./ (b - rho)) - (a * (rho .* rho));
p = p ./ (a * b * b / 27);
plot(rho ./ b, p);

%% Find density for vapor and liquid.
search_range = 0.5:0.001:0.65;
area_delta = [];
last_delta = 100;
search_end = false;
for y = search_range
    root_aux = p - y;
    vapor_dens = find_next_root(root_aux, 0, length(root_aux));
    middle_dens = find_next_root(root_aux, vapor_dens, length(root_aux));
    liquid_dens = find_next_root(root_aux, middle_dens, length(root_aux));
    if (~search_end)
        fprintf('Searching %f, Vid = %d, Mid = %d, Lid = %d, Total = %d\n', y, ...
            vapor_dens, middle_dens, liquid_dens, length(root_aux));
    end

    upper_area = trapz(rho(vapor_dens:middle_dens), root_aux(vapor_dens:middle_dens));
    lower_area = - trapz(rho(middle_dens:liquid_dens), root_aux(middle_dens:liquid_dens));
    this_delta = upper_area - lower_area;
    if (this_delta * last_delta <= 0)
        fprintf('Search End: y = %f, Vapor dens = %f, Liquid dens = %f\n', ...
            y, rho(vapor_dens), rho(liquid_dens));
        search_end = true;
    end
    last_delta = this_delta;
    area_delta = [area_delta; upper_area - lower_area];
end
figure;
plot(search_range, area_delta);

theta
