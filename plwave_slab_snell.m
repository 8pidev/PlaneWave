% Simulation of a plane wave hitting a planar interface using snell's law

lambda = 500; % in nm
n0 = 1.33; % index of medium after objective (immersion oil)
n1 = 1.33; % index of covereslip
n2 = 1.33; % index of medium
NA = 0.8; % Numerical apperture of objective (assumed for n0)
L = 170000; % slab width in nm 170um

NAeff = NA/n0; % Effective NA
WD =  2; % in mm
sinteta = NAeff;
DX = 0; % in nm
DY = tan(asin(NAeff))*WD*1E6; % DX in nm

Nang = 2001;
theta_inc = linspace(-asind(NAeff),asind(NAeff),Nang);
%theta_inc = -70;
Etotal = 0;

for iangle = 1:length(theta_inc)
    
theta = theta_inc(iangle); % in degrees 

c = 3E+8; %c = 3E+8 m/s x 10+9 nm/m x 10-9 ns/s = 3E+8 nm/ns;
k0 = 2*pi*n0/lambda;
w = k0*c/n0;

x = [-3000:50:3000]; % in nm
%x = [-20000:200:20000]; % in nm
y = -WD*1E6 + [-300000:2000:300000]; % in nm
ywater = -2290000:50:-2270000;
yglass = -2310000:50:-2290000;
ywater08 = -1515000:50:-1495000;

y = ywater08;
watermax = [184,61];
watermax08 = [184,61];
glassmax = [160,61];

[X,Y] = meshgrid(x,y);

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

%DX = 0; DY = 0;

% Prepare the new file.
%folder = 'C:\Users\jripoll\OneDrive - Universidad Carlos III de Madrid\UC3M\Powerpoint\';
%vidObj = VideoWriter([folder 'plwave_inf.avi']);
%vidObj.FrameRate = 30;
%open(vidObj);
    
%t = 2*pi/w*(0:0.05:5);
t = 0;
for it=1:length(t)
% %    E = (t01.*t12.*(Y<=-L) + 1.*(Y>-L)).*exp(1i*(k0_x*X + k0_y*Y).*(Y>0) ...
%     E = exp(1i*(k0_x*X + k0_y*Y).*(Y>0) ...
%         + 1i*(k1_x*X + k1_y*Y).*(Y<=0).*(Y>-L) ...
%         + 1i*(k2_x*X + k2_y*Y + (k0_x - k1_x)*L + (k0_y - k1_y)*L).*(Y<=-L) ...
%         ).*exp(1i*w*t(it)).*exp(1i*(k0_x*DX + k0_y*DY));
% %        ).*exp(1i*w*t(it))./(1 + r10.*r12.*exp(-1i*k1_y*L).*(Y<=-L));
% %    E = (t01.*t12.*(Y<=-L) + 1.*(Y>-L)).*exp(1i*(k0_x*X + k0_y*Y).*(Y>0) ...
    E = t01.*t12.*exp(1i*(k2_x*X + k2_y*Y + (k0_x - k1_x)*L + (k0_y - k1_y)*L)) ...
        .*exp(1i*w*t(it)).*exp(1i*(k0_x*DX + k0_y*DY))./(1 + r10.*r12.*exp(-1i*k1_y*L));
    Etotal = Etotal + E;

%         figure(1)
%     hold off
%     subplot(1,2,1)
%     imagesc(x*1E-3,y*1E-3,real(E));
%     xlabel('X (\mum)');
%     ylabel('Y (\mum)');
%     title(['Wavelength in vacuum ' num2str(lambda) 'nm']); 
%     axis xy tight
%     hold on
% %     if n0~=n1
% %         h_hor = line([min(x) max(x)]*1E-3,[0 0],'color','black','LineWidth',0.1,'LineStyle','-');
% %         h_ver = line([0 0],[min(y) max(y)]*1E-3,'color','black','LineWidth',0.1,'LineStyle',':');
% %         h_slab = line([min(x) max(x)]*1E-3,[-L -L]*1E-3,'color','black','LineWidth',0.1,'LineStyle','-');
% %     end
%     subplot(1,2,2)
%     imagesc(x*1E-3,y*1E-3,abs(Etotal).^2);
%     xlabel('X (\mum)');
%     ylabel('Y (\mum)');
%     title(['Wavelength in vacuum ' num2str(lambda) 'nm']); 
%     axis xy tight
%     hold on
% %     if n0~=n1
% %         h_hor = line([min(x) max(x)]*1E-3,[0 0],'color','black','LineWidth',0.1,'LineStyle','-');
% %         h_ver = line([0 0],[min(y) max(y)]*1E-3,'color','black','LineWidth',0.1,'LineStyle',':');
% %         h_slab = line([min(x) max(x)]*1E-3,[-L -L]*1E-3,'color','black','LineWidth',0.1,'LineStyle','-');
% %     end
%     pause(0.1);
% 

%    f = getframe(gcf);
%    writeVideo(vidObj,f);
end



end
    figure
    imagesc(x*1E-3,y*1E-3,abs(Etotal).^2);
    xlabel('X (\mum)');
    ylabel('Z (\mum)');
    title(['Wavelength in vacuum ' num2str(lambda) 'nm']); 
    axis xy tight

% Close the file.
%close(vidObj);