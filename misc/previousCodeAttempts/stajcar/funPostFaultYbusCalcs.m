function [Ybus,Branch] = funPostFaultYbusCalcs(Bus,Branch,Gen,Ybus,Switch)

switch Switch.Type
            
    case 1 % Bus fault.
        Ybus.Y11(Ybus.N,Ybus.N) = Ybus.PostSwitchNN;
        [Ybus] = funYredYrecovCalc(Ybus);

    case 2 % Open branch.

        for jj = 1:length(Branch)
            if (Switch.BusNumNear == Branch(jj).FmBus) && ...
               (Switch.BusNumFar == Branch(jj).ToBus)
                OpenBranchNum = jj;
            elseif (Switch.BusNumNear == Branch(jj).ToBus) && ...
                   (Switch.BusNumFar == Branch(jj).FmBus)
                OpenBranchNum = jj;
            end
        end
        clear jj

        kk = 1;
        for jj = 1:length(Branch)
            if jj ~= OpenBranchNum
                NewBranch(kk) = Branch(jj); 
                kk = kk + 1;
            end
        end
        Branch = NewBranch;
        clear jj kk NewBranch

        clear Ybus
        [Ybus] = funYbusCalc(Bus,Branch,Gen);
        [Ybus] = funBuildYredYrecov(Ybus);
end

end