%% updated class with single get data function
classdef TestClass2 < handle
    properties
        dataNdx = 0;
    end
    methods
%         function obj = TestClass2()
%             % all initializations, calls to base class, etc. here,
%         end
        
        function r = getData(obj, varargin)
            %% getData(k) returns kth data point, if k not provided return all datas
            if obj.dataNdx ~= 0
                global g
                fprintf('num of arguments: %d\n', nargin)
                disp(varargin)
                if nargin == 1
                    disp('Returning:')
                    disp(g.data.list(obj.dataNdx,:))
                    r = g.data.list(obj.dataNdx,:);
                elseif nargin == 2
                    disp('Returning:')
                    disp(g.data.list(obj.dataNdx,varargin{1}))
                    r = g.data.list(obj.dataNdx,varargin{1});
                end
                
            else
                disp('dataNdx not defined')
            end
        end% end getData
        
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

