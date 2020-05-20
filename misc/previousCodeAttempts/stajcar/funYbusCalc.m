function [Ybus] = funYbusCalc(Bus,Branch,Gen)

% Matt Stajcar

% Build power flow Ybus as well as the modified Ybusses that are used for
% solving XXX. The modified Ybusses are described in Power Systems Analysis
% and Design, 5th, by Glover, in chapter 11.5.

% Create power flow Y bus using Branch characteristics only. This 
% calculation is kept separtate from the modified Y bus calculations in 
% order to use this code for a new power flow solver in the future.

Ybus.Y = zeros(length(Bus));
Ybus.Mag = zeros(length(Bus));
Ybus.AngR = zeros(length(Bus)); 

% Calculate normal Y bus.
for row = 1:length(Ybus.Y) 
    for col = row:length(Ybus.Y) 
        % Loop for all rows and columns of upper triangular matrix. This
        % takes advantage of symmetry of Ybus about the diagonal.

        if row == col % Diagonal Ybus element.
            for ii = 1:length(Branch) % Check all branches.
                if (Branch(ii).ToBus == col) || (Branch(ii).FmBus == col)% The 
                    % Branch must be connected to bus in question.
                    Ybus.Y(row,col) = Ybus.Y(row,col) + Branch(ii).Y + Branch(ii).Bpu/2;
                end
            end
        elseif row ~= col % Off diagonal Ybus element.
            for ii = 1:length(Branch) % Check all branches.
                if (Branch(ii).ToBus == col && Branch(ii).FmBus == row) || (Branch(ii).FmBus == col && Branch(ii).ToBus == row)
                    % The Branch must be connected to bus in question.
                    Ybus.Y(row,col) = Ybus.Y(row,col) - Branch(ii).Y;% - Branch(ii).Bpu/2;
                    Ybus.Y(col,row) = Ybus.Y(row,col);
                end
            end
        end
        % Create Ybus with magnitudes and angles.
        Ybus.Mag(row,col) = abs(Ybus.Y(row,col));
        Ybus.AngR(row,col) = angle(Ybus.Y(row,col));

    end
end
clear row col

% Calculate modified Y buses, Y11, Y12, and Y22, as described on page 614
% of the Power System Analysis and Design, 5th, by Glover.

% Y11 - normal Ybus from power flow with generator admittances and load 
% admittances added into it.
Ybus.Y11 = Ybus.Y;
% Add gen admittances to make Y11.
for ii = 1:length(Gen)
    for jj = 1:length(Bus)
        if Bus(jj).Num == Gen(ii).BusNum
            Ybus.Y11(jj,jj) = Ybus.Y11(jj,jj) + Gen(ii).Y;
        end
    end
end
clear ii jj
% Add load admittances to make Y11.
for ii = 1:length(Bus)
    if Bus(ii).Type == 3
        Ybus.Y11(ii,ii) = Ybus.Y11(ii,ii) + Bus(ii).Yload;
    end
end
clear ii

% Y22 - diagonal matrix of generator admittances.
Ybus.Y22 = zeros(length(Gen));
for ii = 1:length(Gen)
    Ybus.Y22(ii,ii) = Gen(ii).Y;
end
clear ii

% Y12 - generator to bus interconnecting admittancesYbu. Y12(k,n) = -gen(n).Y 
% for generator n connected to bus k. If generator n is not connected to
% bus k then Y12(k,n) = 0;
Ybus.Y12 = zeros(length(Bus),length(Gen));
for ii = 1:length(Gen)
    for jj = 1:length(Bus)
        if Bus(jj).Num == Gen(ii).BusNum
            Ybus.Y12(jj,ii) = Ybus.Y12(jj,ii) - Gen(ii).Y;
        end
    end    
end
clear ii jj

end