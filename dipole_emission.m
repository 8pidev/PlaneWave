% from http://pcwww.liv.ac.uk/~awolski/Teaching/Liverpool/PHYS370/AdvancedElectromagnetism-Part8.pdf

lambda = 500;
n0 = 1.0;
n1 = 1.0;
theta = -90; % in degrees 

c = 3E+8; %c = 3E+8 m/s x 10+9 nm/m x 10-9 ns/s = 3E+8 nm/ns;
k0 = 2*pi*n0/lambda;
w = k0*c/n0;
p0 = 1; % normalizing by dipole moment and eps0
eps0 = 1;
x = [-2001:2.5:2000];
y = [-2001:2.5:2000];
%x = [-501:1.5:500];
%y = [-501:1.5:500];
[X,Y] = meshgrid(x,y);

r = sqrt(X.^2 + Y.^2);
theta = unwrap(angle(X + 1i*Y));

t = 2*pi/w*(0:0.05:5);
for it=1:length(t)
    coswt = cos(w*t(it) - k0*r);
    sinwt = sin(w*t(it) - k0*r);
    Er = 1/(4*pi*eps0)*2*k0*p0.*cos(theta)./r.^2.*(sinwt - coswt./(k0.*r)); % real part in far field
    Ephi = -1/(4*pi*eps0)*k0*k0*p0.*sin(theta)./r.*((1-1./(k0.*r).^2).*coswt + sinwt./(k0*r));
    Hphi = -1*c/(4*pi*eps0)*k0*k0*p0.*sin(theta)./r.*(coswt + sinwt./(k0*r));
    Sr = Ephi.*Hphi;
    Sphi = -Er*Hphi;
    Sr(~isfinite(Sr)) = 0;
    Sr(Sr>mean(Sr(:))) = mean(Sr(:));
    Er(~isfinite(Er)) = 0;   
    imagesc(Er.*(abs(r)>100).*r);
    caxis([-0.00001 0.00001])
    pause(0.1)
end
