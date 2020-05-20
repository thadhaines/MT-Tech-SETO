% Matt Stajcar

% This function converts bus, branch, and gen parameters from type
% struct to type array. This enables the user to take advantage of PSTs
% load flow solver.

function [Bus] = Array2Struct(bus_sol,Bus)
    for ii = 1:length(Bus)
        Bus(ii).Vmag = bus_sol(ii,2);
        Bus(ii).AngD = bus_sol(ii,3);
        Bus(ii).Pgen = bus_sol(ii,4);
        Bus(ii).Qgen = bus_sol(ii,5);
        Bus(ii).Pload = bus_sol(ii,6);
        Bus(ii).Qload = bus_sol(ii,7);
    end
    clear ii
end