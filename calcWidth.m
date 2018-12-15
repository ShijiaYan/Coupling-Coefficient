b = 0.2637;
k0 = 2 * pi / 1.55e-6;
nf = 2.397872;
ns = 1.5;
W = zeros(1,4);
disp(' ');
for i = 0:1:3
    V=  (i * pi + 2 * atan(sqrt(b./(1-b)))) ./ sqrt(1-b);
    W(i+1) = V / k0 / (nf^2 - ns^2)^0.5 * 1e6;
    Wprint = sprintf('For i = %d, V = %.3f, W = %.3f microns',i,V,W(i+1));
    disp(Wprint);
    disp(' ');
end