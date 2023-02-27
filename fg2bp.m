% Advanced Orbital Mechanics Assignment 1 Problem 1
% Raman Singh
% Universal Variable two-body orbit propagator Part 2
% Curtis Algorithm 3.4
% All computations are done in metric units

function [r,v] = fg2bp(r_vec0,v_vec0,delta_t,mu)
    % inputs --> initial states, time difference, gravitational parameter

    % magnitudes of the initial state variables
    r0 = norm(r_vec0);
    v0 = norm(v_vec0);

    % radial component of the velocity
    vr0 = dot(r_vec0,v_vec0)/r0;

    alpha = (2/r0) - ((v0^2)/mu);
    
%     if alpha > 0
%         fprintf("Elliptical trajectory")
%     elseif alpha == 0
%         fprintf("Parabolic trajectory")
%     else
%         fprintf("Hyperbolic trajectory")
%     end

    % determining the universal anomaly
    [chi_f] = kepunipro(mu,delta_t,r0,vr0,alpha);

    z = alpha*(chi_f^2);

    % Stumpff functions
    if z > 0
        S = (sqrt(z) - sin(sqrt(z)))/(z^1.5);
        C = (1 - cos(sqrt(z)))/z;
    elseif z == 0
        S = 1/6;
        C = 1/2;
    else
        S = (- sqrt(-z) + sinh(sqrt(-z)))/((-z)^1.5);
        C = (1 - cosh(sqrt(-z)))/z;
    end

    % f & g functions
    f = 1 - ((chi_f^2)/r0)*C;
    g = delta_t - ((1/sqrt(mu))*(chi_f^3)*S);

    r = f*r_vec0 + g*v_vec0;
    
    fdot = ((sqrt(mu))/(norm(r)*r0))*((z*chi_f*S) -chi_f);
    gdot = 1 - ((chi_f^2)/norm(r))*C;

    v = fdot*r_vec0 + gdot*v_vec0;

end