% ldfwdemo.m
% m file to demo loadflow solutions
% a loadflow report is placed in lf_rep.txt if requested

clear
clear global
% global variables
pst_var 
global gen_chg_idx

disp('loadflow demo program')
% load input data from m.file

[dfile,pathname]=uigetfile('d*.m','Select Data File');
if pathname == 0
  error(' you must select a valid data file')
else
  lfile =length(dfile);
  % strip off .m and convert to lower case
  dfile = dfile(1:lfile-2);
  eval(dfile);
end
% check for valid dynamic data file
if isempty(bus)
   error(' the selected file is not a valid data file')
end
answer=input('Do you need a load-flow solution report?  [y/n]n >>     ','s');
if ~isempty(answer)
   if answer == 'Y';answer = 'y';end
   if answer == 'y'; diary lf_rep.txt; end
end
nbus = length(bus(:,1));
gen_chg_idx = ones(nbus,1);
[bus_sol,line_sol,line_flow] = loadflow(bus,line,1e-9,10, ...
                               1.0,answer,2);
if ~isempty(answer)
   if answer=='y'
      diary off
      disp('A load-flow study report is in lf_rep.txt')
   end
end


flag = 0;
while(flag == 0)
  disp('You can examine the system data')
  disp('Type 1 to see initial bus data')
  disp('     2 to see modified line data')
  disp('     3 to see solved load flow bus solution')
  disp('     4 to see line flow')
  disp('     5 to see bus voltage magnitude profile')
  disp('     6 to see bus voltage phase profile')
  disp('     0 to quit')
  sel = input('enter selection >> ');
  if isempty(sel);sel=0;end
  if sel == 1
      disp('                                   Initial Bus Data')
      disp('                                      GENERATION             LOAD')
      disp('       BUS     VOLTS     ANGLE      REAL  REACTIVE      REAL  REACTIVE ')
      disp(bus(:,1:7))
      disp('paused: press any key to continue')
      pause
    elseif sel == 2
      disp('                               Modified Line  Data')
      disp(line_sol)
      disp('paused: press any key to continue')
      pause
    elseif sel == 3
      disp('                                    Solved Bus Data')
      disp('                                      GENERATION             LOAD')
      disp('       BUS     VOLTS     ANGLE      REAL  REACTIVE      REAL  REACTIVE ')
      disp(bus_sol(:,1:7))
      % Custom loadflow output to compare to other sources...
%       fprintf('\t \t \t Solved Bus Data\t \t \t \n')
%       fprintf('\t \t \t \t \t Generation\t Load\n')
%       fprintf(' Bus\t Voltage\t Angle\t Reative\t Real\t Reactive\n')
%       for lNum=1:length(bus_sol)
%         fprintf('%d \t %2.4f \t %2.2f \t %4.2f \t %4.2f \t %4.2f \t %4.2f \n', ...
%             bus_sol(lNum,1),bus_sol(lNum,2),bus_sol(lNum,3),bus_sol(lNum,4)*100, ...
%             bus_sol(lNum,5)*100,bus_sol(lNum,6)*100,bus_sol(lNum,7)*100)
%       end
      disp('paused: press any key to continue')
      pause
    elseif sel == 4
      disp('                  solved Line Flows')
      disp('                      LINE FLOWS                     ')
      disp('      LINE  FROM BUS    TO BUS      REAL  REACTIVE   ')
      disp(line_flow)
      disp('paused: press any key to continue')
      pause
    elseif sel == 5
      bar(bus_sol(:,2))
      title('bus voltage magnitude profile')
      xlabel('internal bus number')
      ylabel('voltage in pu')
      disp('paused: press any key to continue')
      pause
    elseif sel == 6
      bar(bus_sol(:,3))
      title('bus voltage phase profile')
      xlabel('internal bus number')
      ylabel('phase in degrees')
      disp('paused: press any key to continue')
      pause
    elseif sel == 0
      flag = 1;
    else
      error('invalid selection')
  end
end
