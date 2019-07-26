% Simulation of a plane wave hitting a planar interface using snell's law
function E = fun_plwave_slab(theta,lambda,n0,n1,n2,L,DX,DY,X,Y)

c = 3E+8; %c = 3E+8 m/s x 10+9 nm/m x 10-9 ns/s = 3E+8 nm/ns;
k0 = 2*pi*n0/lambda;
k0_x = k0*sind(theta);
k0_y = k0*cosd(theta);

k1 = 2*pi*n1/lambda; % The wavelength in medium1 is lambda/n1
k2 = 2*pi*n2/lambda; % The wavelength in medium1 is lambda/n1
theta1 = asind(k0/k1*sind(theta));
theta2 = asind(k1/k2*sind(theta1));

k1_x = k1*sind(theta1);
k1_y = k1*cosd(theta1);
k2_x = k2*sind(theta2);
k2_y = k2*cosd(theta2);

r10 = (k1_y - k0_y)./(k1_y + k0_y);
r12 = (k1_y - k2_y)./(k1_y + k2_y);
t01 = 2*k0_y./(k1_y + k0_y);
t12 = 2*k1_y./(k1_y + k0_y);

%    E = (t01.*t12.*(Y<=-L) + 1.*(Y>-L)).*exp(1i*(k0_x*X + k0_y*Y).*(Y>0) ...
    E = exp(1i*(X*k2_x + Y*k2_y + ones(size(Y))*(k0_x - k1_x)*L + ones(size(Y))*(k0_y - k1_y)*L) ...
        ).*exp(1i*ones(size(Y))*(k0_x*DX + k0_y*DY));
%        )./(1 + r10.*r12.*exp(-1i*k1_y*L).*(Y<=-L));


return