% Advanced Orbital Mechanics Assignment 1 Problem 2
% Raman Singh
% Lambert Solver
% Curtis Algorithm 5.2
% All computations are done in metric units

function [v1_vec,v2_vec] = lambert_solver(r1_vec,r2_vec,delta_t,mu,orbit_type)
    % inputs --> position vectors at 2 points in space, time difference,
    % gravitational parameter, choice of prograde/retrograde orbit

    % magnitude of the position vectors
    r1 = norm(r1_vec);
    r2 = norm(r2_vec);

    cr = cross(r1_vec,r2_vec);

    % prograde orbit
    switch orbit_type
        % prograde orbit
        case 1
            if cr(3) >= 0
                del_theta = acos(dot(r1_vec,r2_vec)/(r1*r2));
            elseif cr(3) < 0
                del_theta = 2*pi - acos(dot(r1_vec,r2_vec)/(r1*r2));
            end

        % retrograde orbit
        case -1
            if cr(3) < 0
                del_theta = acos(dot(r1_vec,r2_vec)/(r1*r2));
            elseif cr(3) >= 0
                del_theta = 2*pi - acos(dot(r1_vec,r2_vec)/(r1*r2));
            end
    end

    A = sin(del_theta)*sqrt((r1*r2)/(1 - cos(del_theta)));
    
    % tolerance
    tol = 1e-12;

    diff = 1;
    i = 1;
    z(1) = 0.1;

    while diff > tol
        
        % Stumpff function
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

        y(i) = r1 + r2 + A*((z(i)*S(i) - 1)/sqrt(C(i)));
        y0 = r1 + r2 - A*sqrt(2);

        F(i) = ((y(i)/C(i))^1.5)*S(i) + A*sqrt(y(i)) - sqrt(mu)*delta_t;
        
        if z(i) == 0
            Fdot(i) = (sqrt(2)/40)*y0^1.5 + (A/8)*(sqrt(y0) + A*sqrt(1/(2*y0)));
        else
            Fdot(i) = (((y(i)/C(i))^1.5)*((1/(2*z(i)))*(C(i) - 1.5*(S(i)/C(i))) + 0.75*((S(i)^2)/C(i)))) + ...
                ((A/8)*((3*(S(i)/C(i))*sqrt(y(i))) + A*sqrt(C(i)/y(i))));
        end

        z(i+1) = z(i) - F(i)/Fdot(i);

        diff = z(i+1) - z(i);

        i = i + 1;

    end

    zf = z(end);
    
    % Stumpff functions
    if zf > 0
        Sf = (sqrt(zf) - sin(sqrt(zf)))/(zf^1.5);
        Cf = (1 - cos(sqrt(zf)))/zf;
    elseif zf == 0
        Sf = 1/6;
        Cf = 1/2;
    else
        Sf = (- sqrt(-zf) + sinh(sqrt(-zf)))/((-zf)^1.5);
        Cf = (1 - cosh(sqrt(-zf)))/zf;
    end

    yf = r1 + r2 + A*((zf*Sf - 1)/sqrt(Cf));

    % Lagrange f & g functions
    f = 1 - yf/r1;
    g = A*sqrt(yf/mu);
%     fdot = (sqrt(mu)/(r1*r2))*sqrt(yf/Cf)*(zf*Sf-1);
    gdot = 1 - yf/r2;

    v1_vec = (1/g)*(r2_vec - f*r1_vec);
    v2_vec = (1/g)*(gdot*r2_vec - r1_vec);

end