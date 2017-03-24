a = 27;
b = 1;
rtm = 7.2;
filename = '../1/rho0';
field = load(filename);
pressure = rtm * b * field ./ (b - field) - a * field .* field;
laplace = zeros(size(field));
for i = 2:(size(field, 1) - 1)
    for j = 2:(size(field, 2) - 1)
        laplace(i,j) = field(i - 1, j) + field(i + 1, j) + field(i, j - 1) + field(i, j + 1) - 4 * field(i, j);
    end
end
figure;
pcolor(laplace);
colorbar;

figure;
pcolor(pressure);
colorbar;

figure;
pcolor(field);
colorbar;