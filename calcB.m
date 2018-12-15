% This program calculates b.

% Initialization
V = 12.24;
b = zeros(1,4);
neff = zeros(1,4);
Neffh = 2.397872;
syms x



% loop for every mode
for i = 1:1:4
    eqn = ((i-1)*pi + atan(sqrt(x/(1-x))) + atan(sqrt(x/(1-x))))/sqrt(1-x) - V == 0;
    b(i) = vpasolve(eqn , x);
    neff(i) = sqrt(b(i) * (Neffh^2-1.5^2) + 1.5^2);
end
    X1 = sprintf('b = %.4f, %.4f, %.4f, %.4f',b);
    X2 = sprintf('Neff = %.4f, %.4f, %.4f, %.4f',neff);
    disp(" ");
    disp(X1);
    disp(" ");
    disp(X2);
    disp(" ");