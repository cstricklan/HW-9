%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Diffraction Through a Grating
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize MATLAB
close all; clc;
clear all; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Problem
%   
%   
%   
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% UNITS
meters = 1;
decimeters = 1e-1 * meters;
centimeters = 1e-2 * meters;
millimeters = 1e-3 * meters;
inches = 2.54 * centimeters;
feet = 12 * inches;
seconds = 1;
hertz = 1/seconds;
kilohertz = 1e3 * hertz;
megahertz = 1e6 * hertz;
gigahertz = 1e9 * hertz;

% Constants
c0 = 299792458;

% Frequency we want to transmit
f0 = 10 * gigahertz;

NPML = [0 0 20 20];

% Bragg Grating Materials
er1 = 9;

d = 0.75 * centimeters; % the Binary stand height
PeriodWidth = 1.5 * centimeters;  %Width between to Binary stands
dwidth = .5 * PeriodWidth;

%Calculate the Length of our layers.


dc.x = PeriodWidth;
dc.y = d;

Size.x = PeriodWidth;
Size.y = d * 3;

rNx = ceil(Size.x/(centimeters*decimeters)); %This Nz represents real world size
rNy = ceil(Size.y/(centimeters*decimeters));
disp(rNx);
disp(rNy);

% Material Vectors Initialized at Air
rER = ones(rNx,rNy);
rUR = ones(rNx,rNy); 

dx = Size.x/rNx;
dy = Size.y/rNy;


% Fill our Stand
for nx=1:floor(dwidth/dx)
  for ny=1:floor(d/dy)
    rER(nx,ny)=er1;
  end
end


% Fill our body
for nx=1:rNx
  for ny = floor(d/dy)+1:rNy
    rER(nx,ny) = er1;
  end
end


% Add our Materials to the model

% Frequency

freq_start = 0;
freq_end = 15 * gigahertz;

NFREQ = freq_end / (0.5*gigahertz); %Frequencies every 100nm
FREQ = linspace(freq_start, freq_end, NFREQ); %FREQ List

Buffer.x.value = 0;
Buffer.x.e = [-1 -1];
Buffer.x.u = [1 1];

Buffer.y.value = -1;
Buffer.y.e = [1 9];
Buffer.y.u = [1 1];

subplot(121);
imagesc(rER);
% global ERzz;
FDTD2D( dc, Size, rER, rUR, -1, 5e-4, Buffer, NPML, FREQ, NFREQ, 100, 10*gigahertz, 'Binary Grating');
