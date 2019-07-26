% Simulation of a plane wave hitting a planar interface using snell's law

lambda = 500; % in nm
n0 = 1.33; % index of medium after objective (immersion oil)
n1 = 1.5; % index of covereslip
n2 = 1.33; % index of medium
NA = 1.0; % Numerical apperture of objective (assumed for n0)
L = 170000; % slab width in nm 170um

NAeff = NA/n0; % Effective NA
WD =  2; % in mm
sinteta = NAeff;
DX = 0; % in nm
DY = tan(asin(NAeff))*WD*1E6; % DX in nm

Etotal = integral(@(theta_inc) fun_plwave_slab(theta_inc,lambda,n0,n1,n2,L,DX,DY,X(:),Y(:)),...
    -asind(NAeff),asind(NAeff),'ArrayValue',true);
    figure
    imagesc(x*1E-3,y*1E-3,abs(Etotal).^2);
    xlabel('X (\mum)');
    ylabel('Y (\mum)');
    title(['Wavelength in vacuum ' num2str(lambda) 'nm']); 
    axis xy tight

% Close the file.
%close(vidObj);