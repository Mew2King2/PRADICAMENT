%%
%%% analytical example of Fox et al. (2024)

clear all; clc;

% user inputs
b0 = 10;  % T-mm
LB = 1;  % mm, radius of plasma bubble
%dx = 28;  % um, desired/assumed spatial resolution (um/pix)
delta = 30;  % um, half-width of current sheet
KB = 59;  % T, rigidity factor for proton deflection
I0_const = 100;  % protons/pixel, initial source intensity
num_pts = 1000;  % 2 * L_B / dx;  % number of pts

% unit conversions
m_per_mm = 1e-3;
mm_per_m = 1 / m_per_mm;
m_per_um = 1e-6;
um_per_m = 1 / m_per_um;
% 
b0 = b0 * m_per_mm;  % T-m
LB = LB * m_per_mm;  % m
%dx = dx * m_per_um;  % m
delta = delta * m_per_um;  % m

% derived from user inputs
I0 = I0_const * ones(1, num_pts);
% 
x = linspace(-LB, LB, num_pts);
b = Fox_Harris_sheet(b0, delta, LB, x);

% 
squiggle = 1/KB * b;
xprime = x + squiggle;
dxprime = 1 + b0./(KB * delta * LB * abs(x)) .* ( (LB * abs(x) - x.^2).*sech(x/delta).^2 - delta .* x .* tanh(x/delta) );
I_xprime = I0 ./ abs(dxprime);

% 
var_int_term = cumtrapz(x,I0 ./ I_xprime);
cum_int_term = trapz(x,I0 ./ I_xprime);
C = LB - cum_int_term;  % from B.C. b(x = LB = 1) = 0
b_from_inversion = KB * (var_int_term - x + C);

% plot input b(x)
figure;
plot(x * mm_per_m, b * mm_per_m);
xlabel('x (mm)');
ylabel('b(x) (T-mm)');
title('input analytical line-integrated B-field: b(x) = \int dz B_y');

% plot x -> x' mapping
figure;
plot(x * mm_per_m, xprime * mm_per_m);
xlabel('x (mm)');
ylabel('x'' (mm)');
title('calculated x -> x'' mapping: x'' = x + \xi(x)');

% plot I0(x) and I(x')
figure;
hold on;
plot(x * mm_per_m, I0, 'DisplayName', 'I_0(x)');
plot(x * mm_per_m, I_xprime, 'DisplayName', 'I(x'')');
xlabel('x (mm)');
ylabel('fluence (particles/pixel)');
title('forward model for detector image I(x'') from source fluence I_0(x)');
legend(gca, 'show');
hold off;

% plot b (from inversion) and b (analytical input)
figure;
hold on;
plot(x * mm_per_m, b * mm_per_m, 'linewidth', 4, 'DisplayName', 'input b(x)');
plot(x * mm_per_m, b_from_inversion * mm_per_m, 'linewidth', 2, 'DisplayName', 'b(x) from inversion');
xlabel('x (mm)');
ylabel('b(x) (T-mm)');
title('starting b(x) and ending b(x)');
legend(gca, 'show');
hold off;

function b=Fox_Harris_sheet(b0, delta, L_B, x)
b = b0 * tanh(x/delta) .* (1 - abs(x)/L_B);
end
%% 
%%% 