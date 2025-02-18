%%
%%% re-creation of "IV. Verification" section of Fox et al. 2024

clear all; clc;

% user inputs
b0 = 10;  % T-mm
LB = 1;  % mm, radius of plasma bubble
%dx = 28;  % um, desired/assumed spatial resolution (um/pix)
delta = 30;  % um, half-width of current sheet
KB = 59;  % T, rigidity factor for proton deflection
I0_const = 100;  % protons/pixel, initial source intensity

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
num_pts = 1000;  % 2 * L_B / dx;  % number of pts
x0 = linspace(-LB, LB, num_pts);
b=Fox_Harris_sheet(b0, delta, LB, x0);
I0 = I0_const * ones(num_pts,1);

% inputs:
% x0: spatial coordinate where B defined
% b: deflecting line-integrated B field
% KB = deflection parameter (constant)
% I0: undisturbed proton fluence, defined at points x0
%
% outputs:
% X = points in plasma plane where I is determined
% I = proton fluence
[X,I]=prad_fwd(x0,b,KB,I0);

%
figure;
plot(x0 * mm_per_m, b * mm_per_m);

% 
figure;
plot(X * mm_per_m, I);

function b=Fox_Harris_sheet(b0, delta, L_B, x)
b = b0 * tanh(x/delta) .* (1 - abs(x)/L_B);
end
%% 
%%% 