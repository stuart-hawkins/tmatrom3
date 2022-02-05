% TMATROM3 example
%
% This example:
%
% 1) Computes the T-matrix of a sound-hard spheroid
% 2) Visualises the spheroid
% 3) Computes the T-matrix of the spheroid*
% 4) Uses the T-matrix to compute the scattered field for an incident plane wave
% 5) Plots the differential scattering cross section
% 6) Plots the total scattered field
%
% * the T-matrix is saved to disk and if the example is rerun then the
% saved T-matrix is used.
%
% Stuart C. Hawkins - 2 September 2021

% Copyright 2019-2022 Stuart C. Hawkins
% 	
% This file is part of TMATROM3
% 
% TMATROM3 is free software: you can redistribute it and/or modify	
% it under the terms of the GNU General Public License as published by	
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% TMATROM3 is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with TMATROM3.  If not, see <http://www.gnu.org/licenses/>.


%------------------------------------------------
% print info
%------------------------------------------------

help example_sound_hard_spheroid

%------------------------------------------------
% set main parameters
%------------------------------------------------

% set wavenumber
kwave = 5;

% set aspect ratio of prolate spheroid
ar = 1.1;

% set spheroid centre
x0 = [0;0;0];

% set plane wave direction
d = [0;0;1];

% set solver parameter
m = 30;

%------------------------------------------------
% try to load the T-matrix
%------------------------------------------------

% initialise flag to say we didn't load the T-matrix
loaded_tmatrix = 0;

% check if a T-matrix file exists
if exist('tmatrix_sound_hard_spheroid_example.mat','file')
    
    % load the T-matrix
    T = tmatrix('tmatrix_sound_hard_spheroid_example.mat');
    
    % check the wavenumber matches...
    if abs(T.kwave - kwave)<1e-15
   
        % ...if it matched then we accept the saved T-matrix but beware...
        % the aspect ratio or centre could have been changed so it still
        % might not be the correct T-matrix
        
        % update the flag to show we loaded the T-matrix
        loaded_tmatrix = 1;
        
    end
    
end
        
% display a warning message
if loaded_tmatrix
    
    fprintf('Warning:\n\n');
    fprintf('Loaded T-matrix from tmatrix_sound_hard_spheroid_example.mat.\n');
    fprintf('This may not be the correct T-matrix if you changed any of the\n');
    fprintf('experiment parameters.\n\n')
    fprintf('In this case please delete tmatrix_sound_hard_spheroid_example.mat\n');
    fprintf('and rerun this example.\n');
    
end

%------------------------------------------------
% compute the T-matrix
%------------------------------------------------

% setup solver for sound hard scattering by a spheroid
slvr = solverEllipsoidHard(kwave,[],x0,[1 1 ar],m,2*m);
slvr.setup

% work out the suggested wavefunction expansion order
% Note: the diameter of the prolate spheroid is 2*ar so the
% "radius" is ar
n = suggestedorder(kwave,ar);

% compute the T-matrix if we didn't already load it
if ~loaded_tmatrix

    % setup the solver... this assembles the linear system
    slvr.setup
    
    % compute the T-matrix with expansion origin [0;0;0]
    T = ghtmatrix(n,kwave,slvr,[0;0;0]);

    % save the T-matrix
    T.save('tmatrix_sound_hard_spheroid_example.mat')
    
end

%------------------------------------------------
% visualise the scatterer
%------------------------------------------------

% create a new figure
figure(1)

% visualise the scatterer
slvr.visualise()

% make the figure look nice
axis equal
light 
title('Visualisation of the scatterer')

%------------------------------------------------
% use the T-matrix to solve the scattering problem
%------------------------------------------------

% setup incident wave
inc = plane_wave(d,kwave);

% get the regular wavefunction expansion of the incident wave with 
% expansion origin [0;0;0]
a = regularwavefunctionexpansion(n,[0;0;0],inc);

% compute the radiating wavefunction expansion of the scattered wave
% with expansion origin [0;0;0] using the T-matrix
b = T * a;

%------------------------------------------------
% visualise the differential scattering cross section
%------------------------------------------------

% create a new figure
figure

% generate a vector of scattering angles
theta = linspace(0,pi,1000);

% compute corresponding points on the unit sphere at which we
% compute the differential scattering cross section
p = [sin(theta);zeros(size(theta));cos(theta)];

% get the far field at the points in p
psi = b.evaluateFarField(p);

% plot the differential scattering cross section
plot(theta,abs(psi).^2,'b-')

% make the figure look nice
xlabel('scattering angle')
ylabel('differential scattering cross section')
title('Differential scattering cross section')

%------------------------------------------------
% visualise the total field
%------------------------------------------------

% create a new figure
figure

% create a mesh of points at which to compute the field
u = linspace(-3,3,100);
v = linspace(-3,3,100);
[x,z] = meshgrid(u,v);
y = zeros(size(x));

% stack the mesh points into a 3 x m array
p = [x(:).';y(:).';z(:).'];

% the wavefunction expansion is only valid outside a sphere of radius
% ar... we need to apply a mask
r = sqrt(sum(p.^2,1));
mask = (r>1.05*ar);

% initialise an array to hold the field
psi = zeros(size(x));

% get the field at the masked mesh points
tmp = b.evaluate(p(:,mask)) + inc.evaluate(p(:,mask));

% slot the field values into their array
psi(mask) = tmp;

% make the masked points NaN so that they don't get plotted
psi(~mask) = NaN;

% plot the real part of the field
surf(x,y,z,real(psi))

% make the figure look nice
axis equal
shading interp
title('Total field (real part)')

% add the scatterer
hold on
slvr.visualise();
hold off