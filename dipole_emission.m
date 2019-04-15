% from http://pcwww.liv.ac.uk/~awolski/Teaching/Liverpool/PHYS370/AdvancedElectromagnetism-Part8.pdf

circd = @(radius,deg_ang)  [radius*cosd(deg_ang);  radius*sind(deg_ang)];

lambda = 500;
n0 = 1.0;
n1 = 1.0;
theta = -90; % in degrees 

c = 3E+8; %c = 3E+8 m/s x 10+9 nm/m x 10-9 ns/s = 3E+8 nm/ns;
k0 = 2*pi*n0/lambda;
w = k0*c/n0;
p0 = 1; % normalizing by dipole moment and eps0
eps0 = 1e-3;
x = [-2000:5:2000];
y = [-2000:5:2000];
%x = [-500:1:500];
%y = [-500:1:500];
[X,Y] = meshgrid(x,y);

r0 = sqrt(X.^2 + Y.^2);
theta = unwrap(angle(X + 1i*Y));

fact = logspace(-3,0,10);

for ifact = 1:length(fact)
    r = r0*fact(ifact);
    
    L = 100; Rc = 20;
    t = 2*pi/w*(0:0.05:1);
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
        figure(1)
        hold off
        imagesc(Er.*r.^2);
        caxis([-1/fact(ifact) 1/fact(ifact)])
        hold on
        contour(real(log(Er.*r^2)),'k:')
        ypos = 395 + cos(w*t(it))*L;
        circ = circd(Rc,linspace(0,360,90));
        fill(398+circ(1,:),ypos+circ(2,:),'w');
        text(398+6-Rc,ypos+6-Rc,'-','FontSize',30)
        pause(0.1)
    end
end