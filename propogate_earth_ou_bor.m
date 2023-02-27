% Advanced Orbital Mechanics Assignment 1 Problem 1
% Raman Singh
% Universal Variable two-body orbit propagator Part 3
% All computations are done in metric units

close all; clear; clc;
AU = 149597870.7; % km
day = 60*60*24; % seconds

% initial states for Earth
RiE = [-1.796136509111975e-1,9.667949206859814e-1,-3.668681017942158e-5]*AU; % in AU
ViE = [-1.720038360888334e-2,-3.211186197806460e-3,7.927736735960840e-7]*(AU/day);% in AU/day

% initial states for Oumouamoua
Ri1I = [3.515868886595499e-2,-3.162046390773074,4.493983111703389]*AU; % in AU
Vi1I = [-2.317577766980901e-3,9.843360903693031e-3,-1.541856855538041e-2]*(AU/day); % in AU/day

% initial states for Borisov
Ri2I = [7.249472033259724,14.61063037906177,14.24274452216359]*AU; % in AU
Vi2I = [-8.241709369476881e-3,-1.156219024581502e-2,-1.317135977481448e-2]*(AU/day); % in AU/day

tint = linspace(0,3.5*365*day,3e4);
mu = 1.3271244e11; % gravitational parameter for the Sun

% propagating Earth's trajectory
pro_states_E = zeros(length(tint),6);

for i = 1:length(tint)
    [pro_states_E(i,1:3),pro_states_E(i,4:6)] = fg2bp(RiE,ViE,tint(i),mu);
end

% propagating Oumuoamuoa's trajectory
pro_states_O = zeros(length(tint),6);

for i = 1:length(tint)
    [pro_states_O(i,1:3),pro_states_O(i,4:6)] = fg2bp(Ri1I,Vi1I,tint(i),mu);
end

% propagating Borisov's trajectory
pro_states_B = zeros(length(tint),6);

for i = 1:length(tint)
    [pro_states_B(i,1:3),pro_states_B(i,4:6)] = fg2bp(Ri2I,Vi2I,tint(i),mu);
end

% plotting
figure
hold on

plot3(pro_states_E(:,1)/AU,pro_states_E(:,2)/AU,pro_states_E(:,3)/AU,Color='b',LineWidth=1)
plot3(pro_states_O(:,1)/AU,pro_states_O(:,2)/AU,pro_states_O(:,3)/AU,Color='k',LineWidth=1)
plot3(pro_states_B(:,1)/AU,pro_states_B(:,2)/AU,pro_states_B(:,3)/AU,Color='m',LineWidth=1)
plot3(0,0,0,'r*');

legend([{"Earth"},{"Oumuoamuoa"},{"Borisov"},{"Sun"}])
xlabel('x in AU')
ylabel('y in AU')
zlabel('z in AU')
title('Orbit Plot for Earth, Oumuoamuoa and Borisov')
view([-15 30])