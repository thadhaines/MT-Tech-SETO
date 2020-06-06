%% simple test class with data get and set to global g struct
classdef TestClass < handle
    properties
        dataNdx = 0;
    end
    methods
        function obj = TestClass()
            % all initializations, calls to base class, etc. here,
        end
        
        function r = getAllData(obj)
            if obj.dataNdx ~= 0
                global g
                disp('Returning:')
                disp(g.data.list(obj.dataNdx,:))
                r = g.data.list(obj.dataNdx,:);
            else
                disp('dataNdx not defined')
            end
        end% end getAllData
        
        function r = getkData(obj,k)
            if obj.dataNdx ~= 0
                global g
                disp('Returning:')
                disp(g.data.list(obj.dataNdx,k))
                r = g.data.list(obj.dataNdx,k);
            else
                disp('dataNdx not defined')
            end
        end% end getkData
        
        function setData(obj, k , data)
            if obj.dataNdx ~= 0
                global g
                disp('Data Changed from:')
                disp(g.data.list(obj.dataNdx,:))
                
                g.data.list(obj.dataNdx,k) = data;
                disp('To:')
                disp(g.data.list(obj.dataNdx,:))
            else
                disp('dataNdx not defined')
            end
        end% end set data
        
    end% end methods
end % end class def

