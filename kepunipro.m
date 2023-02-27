% Advanced Orbital Mechanics Assignment 1 Problem 1
% Raman Singh
% Universal Variable two-body orbit propagator Part 1
% Curtis Algorithm 3.3
% All computations are done in metric units

function [chi_f] = kepunipro(mu,del_t,r_0,v_r0,alpha)
    % inputs ---> gravitational parameter, time difference, initial
    % position, radial component of velocity, inverse of semi-major axis

    % tolerance
    tol = 1e-9;

    % initial guess for chi
    chi_0 = sqrt(mu)*abs(alpha)*del_t;
    chi(1) = chi_0;
    
    for i = 1:10000
    
        z(i) = alpha*(chi(i)^2);
    
        % Stumpff functions
        if z(i) > 0
            S(i) = (sqrt(z(i)) - sin(sqrt(z(i))))/(z(i)^1.5);
            C(i) = (1 - cos(sqrt(z(i))))/z(i);
        elseif z(i) == 0
            S(i) = 1/6;
            C(i) = 1/2;
        else
            S(i) = (- sqrt(-z(i)) + sinh(sqrt(-z(i))))/((-z(i))^1.5);
            C(i) = (1 - cosh(sqrt(-z(i))))/z(i);
        end
    
        f(i) = (((r_0*v_r0)/sqrt(mu))*(chi(i)^2)*C(i)) + ((1 - alpha*r_0)*(chi(i)^3)*S(i)) + r_0*chi(i) - sqrt(mu)*del_t;
        f_p(i) = (((r_0*v_r0)/sqrt(mu))*chi(i)*(1 - alpha*(chi(i)^2)*S(i))) + ((1 - alpha*r_0)*(chi(i)^2)*C(i)) + r_0;

        ratio(i) = f(i)/f_p(i);
    
        chi(i+1) = chi(i) - ratio(i);
        
        if abs(ratio(i)) < tol
            chi_f = chi(i);
            break;
        end
    end
end