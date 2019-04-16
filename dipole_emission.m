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
%x = [-2000:5:2000];
%y = [-2000:5:2000];
x = [-lambda/2:0.5:lambda/2];
y = [-lambda/2:0.5:lambda/2];
[X,Y] = meshgrid(x,y);

r0 = sqrt(X.^2 + Y.^2);
theta = unwrap(angle(X + 1i*Y));

fact = [1 2 5 10];
showpart = true;

% Prepare the new file.
folder = 'C:\Users\jripoll\OneDrive - Universidad Carlos III de Madrid\UC3M\Powerpoint\';


for ifact = 1:length(fact)
    r = r0*fact(ifact);
    vidObj = VideoWriter([folder 'dipole_emm' num2str(fact(ifact)) '.avi']);
    vidObj.FrameRate = 30;
    open(vidObj);
    
    if fact(ifact)>1, showpart = false; end
    
    L = 10; Rc = 10;
    t = 2*pi/w*(0:0.01:4);
    for it=1:length(t)
        coswt = cos(w*t(it) - k0*r);
        sinwt = sin(w*t(it) - k0*r);
        Er = 1/(4*pi*eps0)*2*k0*p0.*cos(theta)./r.^2.*(sinwt - coswt./(k0.*r)); % real part in far field
        Ephi = -1/(4*pi*eps0)*k0*k0*p0.*sin(theta)./r.*((1-1./(k0.*r).^2).*coswt + sinwt./(k0*r));
        Hphi = -1*c/(4*pi*eps0)*k0*k0*p0.*sin(theta)./r.*(coswt + sinwt./(k0*r));
        Sr = Ephi.*Hphi;
        Sphi = -Er*Hphi;
        Sr(~isfinite(Sr)) = 0;
        Sr(Sr>mea+
        +n(Sr(:))) = mean(Sr(:));
        Er(~isfinite(Er)) = 0; 
        figure(1)
        hold off
        imagesc(x*fact(ifact),y*fact(ifact),Er.*r.^2);
        caxis([-2 2])
        axis xy
        hold on
        contour(x*fact(ifact),y*fact(ifact),real(log(Er.*r^2)),'k:')
        axis off
        text(max(x)*fact(ifact)*0.3,max(y)*fact(ifact)*0.7,['Scale ' num2str(fact(ifact)) '\lambda'],'FontSize',20,'Color','w')
            ypos = - sin(w*t(it))*L;
            circ = circd(Rc*fact(ifact),linspace(0,360,90));
            fill(circ(1,:),ypos+circ(2,:),'w');
        if showpart
            text(-5*fact(ifact),ypos+5*fact(ifact),'-','FontSize',30)
        end
        line(lambda/2*[-1 1 1 -1 -1],lambda/2*[-1 -1 1 1 -1],'color','black','LineWidth',0.1,'LineStyle','-')
        pause(0.001)
        f = getframe(gcf);
        writeVideo(vidObj,f);
    end
% Close the file.
close(vidObj);
end
