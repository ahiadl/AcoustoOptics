% Mie theory
% Use SI units
% Wang

% diameter = input('Diameter of sphere (e.g., 579 nm):')*1e-9;
% radius = diameter/2;
% lambda = input('Wavelength (e.g., 400 nm):')*1e-9;

diameter = 1.15e-6;
radius = diameter/2;
lambda = 780e-9;

% n_s = input('Refractive index of sphere (e.g., 1.57):');
% n_b = input('Refractive index of background (e.g., 1.33):');
% w_s = input('Specific weight of sphere(e.g., 1.05 g/cc):')*1e3;
% w_b = input('Specific weight of background(e.g., 1.0 g/cc):')*1e3;

n_s = 2.5249;
n_b = 1.4;
w_s = 4*1e3;
w_b = 1*1e3;

% concentration = input('Concentration by weight (e.g., 0.002):');

concentration = 0.01;

k     = 2*pi*n_b/lambda;
x     = k*radius;
n_rel = n_s/n_b;
y     = n_rel*x;

% Calculate the summations
err = 1e-8;
Qs = 0;
gQs = 0;

for n = 1:100000
    Snx = sqrt(pi*x/2)*besselj(n+0.5,x);
    Sny = sqrt(pi*y/2)*besselj(n+0.5,y);
    Cnx = -sqrt(pi*x/2)*bessely(n+0.5,x);
    Zetax = Snx+1i*Cnx;
    % Calculate the first-order derivatives
    Snx_prime = - (n/x)*Snx+sqrt(pi*x/2)*besselj(n-0.5,x);
    Sny_prime = - (n/y)*Sny+sqrt(pi*y/2)*besselj(n-0.5,y);
    Cnx_prime = - (n/x)*Cnx-sqrt(pi*x/2)*bessely(n-0.5,x);
    Zetax_prime = Snx_prime + 1i*Cnx_prime;
    an_num = Sny_prime*Snx-n_rel*Sny*Snx_prime;
    an_den = Sny_prime*Zetax - n_rel*Sny*Zetax_prime;
    an = an_num/an_den;
    bn_num = n_rel*Sny_prime*Snx-Sny*Snx_prime;
    bn_den = n_rel*Sny_prime*Zetax-Sny*Zetax_prime;
    bn = bn_num/bn_den;
    Qs1 = (2*n+1)*(abs(an)^2+abs(bn)^2);
    Qs = Qs+Qs1;

    if  n > 1
        gQs1 = (n-1)*(n+1)/n*real(an_1*conj(an)+bn_1*conj(bn))...
        +(2*n-1)/((n-1)*n)*real(an_1*conj(bn_1));
        gQs = gQs+gQs1;
    end
    an_1 = an;
    bn_1 = bn;
    if abs(Qs1)<(err*Qs) && abs(gQs1)<(err*gQs)
    break;
    end
end
Qs = (2/x^2)*Qs;
gQs = (4/x^2)*gQs;
g = gQs/Qs;
vol_s = 4*pi/3*radius^3;
N_s = concentration*w_b/(vol_s*w_s);
sigma_s = 05*pi*radius^2;
mu_s = N_s*sigma_s;
mu_s_prime = mu_s*(1-g);
% Output results
{'wavelength(nm)','Qs ( - )',' g (-)','mus (/cm)','mus_prime(/cm)';...
lambda*1e9, Qs,g, mu_s*1e-2, mu_s_prime*1e-2}