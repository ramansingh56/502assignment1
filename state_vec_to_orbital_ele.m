% Advanced Orbital Mechanics Assignment 1 Problem 5
% Raman Singh
% State vectors to Orbital elements
% Curtis Algorithm 4.2
% All computations are done in metric units

close all; clear; clc;

AU = 149597870.7;
day = 60*60*24;
mu = 1.3271244e11;

% initial states for Earth
RiE = [-1.796136509111975e-1,9.667949206859814e-1,-3.668681017942158e-5]*AU; % in AU
ViE = [-1.720038360888334e-2,-3.211186197806460e-3,7.927736735960840e-7]*(AU/day);% in AU/day

object = 2;

switch object

    case 1
        % initial states for Oumouamoua
        RiA = [3.515868886595499e-2,-3.162046390773074,4.493983111703389]*AU; % in AU
        ViA = [-2.317577766980901e-3,9.843360903693031e-3,-1.541856855538041e-2]*(AU/day); % in AU/day
        
    case 2
        % initial states for Borisov
        RiA = [7.249472033259724,14.61063037906177,14.24274452216359]*AU; % in AU
        ViA = [-8.241709369476881e-3,-1.156219024581502e-2,-1.317135977481448e-2]*(AU/day); % in AU/day

end

[a,e,i,omega,w,f] = P5_rv2oe(RiA,ViA,mu)

function [a,e,i,omega,w,f] = P5_rv2oe(r_vec,v_vec,mu)

    % Convert position & velocity to orbital elements

    i_vec = [1 0 0];
    j_vec = [0 1 0];
    k_vec = [0 0 1];
    r = norm(r_vec);
    v = norm(v_vec);
    
    % vis-viva equation for semi-major axis (a)

    a = 1/((2/r)-((v^2)/mu));

    % eccentricity

    e_vec = (((v^2)/mu)-(1/r))*r_vec - ((1/mu)*(dot(r_vec,v_vec))*v_vec);
    e = norm(e_vec);

    % angular momentum

    h_vec = cross(r_vec,v_vec);
    h = norm(h_vec);

    % inclination
    
    i = acos(dot((h_vec/h),k_vec));

    % node vector

    n_vec = cross(k_vec,h_vec);
    n = norm(n_vec);

    % Longitude of the ascending node

    omega = acos((dot(n_vec,i_vec))/n);

    if dot(n_vec,j_vec) < 0
        omega = 2*pi - omega;
    end

    % argument of periapse

    w = acos((dot(n_vec,e_vec))/(n*e));

    if dot(e_vec,k_vec) < 0
        w = 2*pi - w;
    end

    % true anomaly

    f = acos((dot(r_vec,e_vec))/(r*e));

    if dot(r_vec,v_vec) < 0
        f = 2*pi - f;
    end

end