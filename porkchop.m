% Advanced Orbital Mechanics Assignment 1 Problem 3
% Raman Singh
% Pork Chop Plots
% All computations are done in metric units

close all; clear; clc;

AU = 149597870.7;
day = 60*60*24;
mu = 1.3271244e11;

% initial states for Earth
RiE = [-1.796136509111975e-1,9.667949206859814e-1,-3.668681017942158e-5]*AU; % in AU
ViE = [-1.720038360888334e-2,-3.211186197806460e-3,7.927736735960840e-7]*(AU/day); % in AU/day

object = 1;

switch object

    case 1
        % for Oumouamoua
        % initial states
        RiA = [3.515868886595499e-2,-3.162046390773074,4.493983111703389]*AU; % in AU
        ViA = [-2.317577766980901e-3,9.843360903693031e-3,-1.541856855538041e-2]*(AU/day); % in AU/day
        
        % arrival & departure time intervals
        t_dep = linspace(1*day,365*day,365);
        t_arr = linspace(212*day,761*day,550);

    case 2
        % for Borisov
        % initial states
        RiA = [7.249472033259724,14.61063037906177,14.24274452216359]*AU; % in AU
        ViA = [-8.241709369476881e-3,-1.156219024581502e-2,-1.317135977481448e-2]*(AU/day); % in AU/day
        
        % arrival & departure time intervals
        t_dep = linspace(1*day,1308*day,1308);
        t_arr = linspace(882*day,1857*day,976);
end

mission_type = 2;

switch mission_type
    case 1
        % for rendezvous
        for i = 1:length(t_dep)
            [RfE,VfE] = fg2bp(RiE,ViE,t_dep(i),mu);
            for j = 1:length(t_arr)
                if t_dep(i) - t_arr(j) < 0
                    [RfA,VfA] = fg2bp(RiA,ViA,t_arr(j),mu);
                    [vsE,vsA] = lambert_solver(RfE,RfA,t_arr(j)-t_dep(i),mu,1);
                    del_v(j,i) = norm(vsE - VfE) + norm(VfA - vsA);
                else
                    del_v(j,i) = NaN;
                end
            end
        end

    case 2
        % for flyby
        for i = 1:length(t_dep)
            [RfE,VfE] = fg2bp(RiE,ViE,t_dep(i),mu);
            for j = 1:length(t_arr)
                if t_dep(i) - t_arr(j) < 0
                    [RfA,VfA] = fg2bp(RiA,ViA,t_arr(j),mu);
                    [vsE,vsA] = lambert_solver(RfE,RfA,t_arr(j)-t_dep(i),mu,1);
                    del_v(j,i) = norm(vsE - VfE);
                else
                    del_v(j,i) = NaN;
                end
            end
        end
end

switch mission_type
    case 1
        if object == 1
            for i = 1:length(t_dep)
                for j = 1:length(t_arr)
                    if del_v(j,i) > 50
                        del_v(j,i) = NaN;
                    end
                end
            end
        
        elseif object == 2
            for i = 1:length(t_dep)
                for j = 1:length(t_arr)
                    if del_v(j,i) > 60
                        del_v(j,i) = NaN;
                    end
                end
            end
        end

    case 2
         for i = 1:length(t_dep)
            for j = 1:length(t_arr)
                if del_v(j,i) > 20
                    del_v(j,i) = NaN;
                end
            end
         end
end

% to determine the minimum delta v value
min_del_v = min(min(del_v));
[row_pos,col_pos] = find(del_v==min_del_v);

% plotting
figure
surf(t_dep/day,t_arr/day,del_v,'EdgeColor','none');
colormap('jet')
view([0 90])
hold on
colorbar

switch object
    case 1
        dateaxis('x',12,datetime(2017,1,1))
        dateaxis('y',12,datetime(2017,8,1))
        xlabel('Departure time');
        ylabel('Arrival time');
        if mission_type == 1
            title('Porkchop Plots (\Deltav variation) for rendezvous with Oumouamoua');
        elseif mission_type == 2
            title('Porkchop Plots (\Deltav variation) for fly-by - Oumouamoua');
        end
    case 2
        dateaxis('x',12,datetime(2017,1,1))
        dateaxis('y',12,datetime(2019,6,1))
        xlabel('Departure time');
        ylabel('Arrival time');
        if mission_type == 1
            title('Porkchop Plots (\Deltav variation) for rendezvous with Borisov');
        elseif mission_type == 2
            title('Porkchop Plots (\Deltav variation) for fly-by - Borisov');
        end
end