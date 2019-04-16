% Simulation of a plane wave hitting a planar interface using snell's law

circd = @(radius,deg_ang,xo,yo)  [radius*cosd(deg_ang)+xo;  radius*sind(deg_ang)+yo];

lambda = 500;
n0 = 1.0;
n1 = 1.333;
theta_ref = 0;

c = 3E+8; %c = 3E+8 m/s x 10+9 nm/m x 10-9 ns/s = 3E+8 nm/ns;
k0 = 2*pi*n0/lambda;
w = k0*c/n0;

x = [-2000:2:2000];
y = [-2000:2:2000];
[X,Y] = meshgrid(x,y);

k0_x = k0*sind(theta_ref);
k0_y = k0*cosd(theta_ref);

% The lens arc:
rcenter = [0 3000];
XYlens = circd(3000,linspace(180,360,101),rcenter(1),rcenter(2));  

theta_inc = 90-180/pi*angle((XYlens(1,:)-rcenter(1)) + 1i*(XYlens(2,:)-rcenter(2)));

Rott0 = zeros(size(X));
for ii=10:10:length(theta_inc)
    Rott = (X-XYlens(1,ii))*sind(theta_inc(ii))+(Y-XYlens(2,ii))*cosd(theta_inc(ii));
    Rott0 = Rott0 + (Rott>0);
    figure(1)
    hold off
    imagesc(x,y,Rott0);
    axis xy
    hold on
    plot(XYlens(1,:),XYlens(2,:),'w.')
    pause(0.1)
    
    xn0 = min(x);
    xn1 = max(x);
    yn0 = XYlens(2,ii) + (xn0-XYlens(1,ii))*tand(-theta_inc(ii));
    yn1 = XYlens(2,ii) + (xn1-XYlens(1,ii))*tand(-theta_inc(ii));
    line([xn0 xn1],[yn0 yn1],'color','k');
        
theta = -(180-theta_inc(ii));
k1 = 2*pi*n1/lambda; % The wavelength in medium1 is lambda/n1
theta1 = asind(k0/k1*sind(theta));
k1_x = k1*sind(theta1);
k1_y = k1*cosd(theta1);

arci = circd(500,linspace(90-theta,90,10),XYlens(1,ii),XYlens(2,ii));  


% Prepare the new file.
%folder = 'C:\Users\jripoll\OneDrive - Universidad Carlos III de Madrid\UC3M\Powerpoint\';
%vidObj = VideoWriter([folder 'plwave_inf.avi']);
%vidObj.FrameRate = 30;
%open(vidObj);
    
t = 2*pi/w*(0:0.05:2);
for it=1:length(t)
    new_angle = theta_inc(ii) - theta1;
    k1_xn = k1*sind(new_angle);
    k1_yn = k1*cosd(new_angle);
    
    deph1 = k1_xn*XYlens(1,ii) + k1_yn*XYlens(2,ii);
    deph0 = k0_x*XYlens(1,ii) + k0_y*XYlens(2,ii);
    E = exp(1i*(k0_x*X + k0_y*Y - deph0).*(Rott<=0) - 1i*(k1_xn*X + k1_yn*Y - deph1).*(Rott>0)).*exp(1i*w*t(it));
    
    xn0 = min(x);
    xn1 = max(x);
    yn0 = XYlens(2,ii) + (xn0-XYlens(1,ii))*tand(-theta_inc(ii));
    yn1 = XYlens(2,ii) + (xn1-XYlens(1,ii))*tand(-theta_inc(ii));
    yp0 = XYlens(2,ii) + (xn0-XYlens(1,ii))*tand(90-theta_inc(ii));
    yp1 = XYlens(2,ii) + (xn1-XYlens(1,ii))*tand(90-theta_inc(ii));
   
    figure(2)
    hold off
    imagesc(x,y,real(E));
    xlabel('X (\mum)');
    ylabel('Y (\mum)');
    title(['Wavelength in vacuum ' num2str(lambda) 'nm. Rel-angle ' num2str(theta) 'o']); 
    axis xy
    hold on
    if n0~=n1
        plot(arci(1,:),arci(2,:),'w-','LineWidth',1)
        plot(XYlens(1,:),XYlens(2,:),'w-')
        h_hor = line([xn0 xn1],[yn0 yn1],'color','k');
        h_ver = line([xn0 xn1],[yp0 yp1],'color','black','LineWidth',0.1,'LineStyle',':');
        h_inc = line([XYlens(1,ii) XYlens(1,ii)],[XYlens(2,ii) max(y)],'color','w','LineWidth',0.2,'LineStyle','-'); 
    end
    %text(max(x)*0.3,max(y)*0.7,[],'FontSize',20,'Color','w')
       
    pause(0.1);
%    f = getframe(gcf);
%    writeVideo(vidObj,f);
end

end
% Close the file.
%close(vidObj);