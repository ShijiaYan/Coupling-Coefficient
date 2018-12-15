% This code is to calculate coupling coefficient for different waveguide
% settings.
disp(' ');
% Initialization
c = 3e8;
lambda = 1.55e-6;
h = 0.22e-6;
eps = 8.85e-12;
omega = 2 * pi * c / lambda;
k0 = 2 * pi / lambda;
nf = 3;
ns = 1.5;
mu = 4 * pi * 1e-7;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           INPUT             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wsg = 0.1657e-6;

% Wmx = 0.1657e-6;
% nu = 0;
% Wmx = 0.6485e-6;
% nu = 1;
% Wmx = 1.1313e-6;
% nu = 2;
Wmx = 1.6141e-6;
nu = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       CONFIGURATION         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step = 1e4;
Gap = 1e-6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






% Calculate Neffh
V = k0 * h * sqrt(nf^2 - ns^2);
syms bb
eqn1 =  2 * atan(sqrt(bb/(1-bb)))/sqrt(1-bb) - V == 0;
b = vpasolve(eqn1 , bb);
Neffh = sqrt(b * (nf^2-ns^2) + ns^2); % Good

% Calculate Mode for SMWg
V = k0 * Wsg * (Neffh^2 - ns^2)^0.5;
eqn2 =  2 * atan(sqrt(bb/(1-bb)))/sqrt(1-bb) - V == 0;
b = vpasolve(eqn2 , bb);
Neff1 = sqrt(b * (Neffh^2-ns^2) + ns^2); % Good
% Calculate mode constants
beta1 = k0 * Neff1;
kappa1 = k0 * sqrt(Neffh^2 - Neff1^2);
gamma1 = k0 * sqrt(Neff1^2 - ns^2);

% Calculate the Mode for MMwg
V = k0 * Wmx * (Neffh^2 - ns^2)^0.5;
eqn3 =  (pi * nu + 2 * atan(sqrt(bb/(1-bb))))/sqrt(1-bb) - V == 0;
b = vpasolve(eqn3 , bb);
Neff2 = sqrt(b * (Neffh^2-ns^2) + ns^2); %Good
% Calculate mode constants
beta2 = k0 * Neff2;
kappa2 = k0 * sqrt(Neffh^2 - Neff2^2);
gamma2 = k0 * sqrt(Neff2^2 - ns^2);
disp('Calculate Mode Const Complete.');
Mprint = sprintf('Neff1 = %.4f, Neff2 = %.4f',Neff1,Neff2);
disp(Mprint);
disp(' ');

% Normalization for SMWg
P = 0;
xb = 5e-6;
x = linspace(-xb,xb,step);
% Calculate E field
E1 = zeros(1,step);
dx = x(2) - x(1);
for i = 1:1:step
    if x(i) < -Wsg/2
        E1(i) = exp(gamma1 * ( Wsg/2 + x(i) ));
    end
    if x(i) > Wsg/2
        E1(i) = exp(-gamma1 * ( -Wsg/2 + x(i) ));
    end
    if x(i) >= -Wsg/2 && x(i) <= Wsg/2
        E1(i) = cos(kappa1 * x(i)) / cos(kappa1 * Wsg / 2);
    end
    P = P + E1(i)^2 * dx;
end
% Calculate C
C1 = 1 / sqrt(beta1 * P / 2 / omega / mu);
disp('Normalization for SMW complete.');
C1print = sprintf('C1 is equal to %.4e',C1);
disp(C1print);
disp(' ');

% Normalization for MMWg
% Calculate E field
P = 0;
E2 = zeros(1,step);
for i = 1:1:step
    if nu == 0 || nu == 2
        if x(i) < -Wmx/2
            E2(i) = exp(gamma2 * ( Wmx/2 + x(i) ));
        end
        if x(i) > Wmx/2
            E2(i) = exp(-gamma2 * ( -Wmx/2 + x(i) ));
        end
        if x(i) >= -Wmx/2 && x(i) <= Wmx/2
            E2(i) = cos(kappa2 * x(i)) / cos(kappa2 * Wmx / 2);
        end
    else
        if x(i) < - Wmx/2
            E2(i) = exp(gamma2 * ( x(i) + Wmx/2 ));
        end
        if x(i) > Wmx/2
            E2(i) = -exp(-gamma2 * ( x(i) - Wmx/2 ));
        end
        if x(i) >= -Wmx/2 && x(i) <= Wmx/2
            E2(i) = -sin(kappa2 * x(i)) / sin(kappa2 * Wmx / 2);
        end
    end
    P = P + E2(i)^2 * dx;
end
C2 = 1 / sqrt(beta2 * P / 2 / omega / mu);
disp('Normalization for MMW complete.');
C2print = sprintf('C2 is equal to %.4e',C2);
disp(C2print);
disp(' ');

% Calculate coupling efficiency and Coupling Length
xcoup = linspace(-Wmx/2,Wmx/2,step);
dxcoup = xcoup(2) - xcoup(1);
Int = 0;
for i = 1:1:step
    % Calculate E field
    E1(i) = C1 * exp(-gamma1 * ( xcoup(i) + Gap + Wmx/2 ));
    if nu == 0 || nu == 2
        E2(i) = C2 * cos(kappa2 * xcoup(i)) / cos(kappa2 * Wmx / 2);
    else
        E2(i) = -C2 * sin(kappa2 * xcoup(i)) / sin(kappa2 * Wmx / 2);
    end
    % Calculate the integral
    Int = Int + E1(i) * E2(i) * dxcoup;
end
K = eps * omega * (Neffh^2 - ns^2) * Int / 4;
Lpi = pi / 2 / K;
Lprint = sprintf('Lpi = %.3f microns',Lpi*1e6);
disp(Lprint);
disp(' ');