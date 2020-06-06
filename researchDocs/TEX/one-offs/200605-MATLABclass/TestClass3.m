%% updated class with single get data function, no debug messages...

% assumes a global g has been declared with data in a g.data.list array
classdef TestClass3 < handle %inherit from handle class
    properties
        dataNdx = 0; % placeholder
    end
    methods
        function r = getData(obj, varargin)
            %% getData(k) returns kth data point, if k not provided return all datas
            if obj.dataNdx ~= 0
                global g
                if nargin == 1
                    r = g.data.list(obj.dataNdx,:);
                elseif nargin > 1
                    r = g.data.list(obj.dataNdx,varargin{1});
                end
            else
                disp('dataNdx not defined')
            end
        end% end getData
        
        function setData(obj, k , data)
            %% setData(k, newData) sets kth data to new data
            if obj.dataNdx ~= 0
                global g
                g.data.list(obj.dataNdx,k) = data;
            else
                disp('dataNdx not defined')
            end
        end% end setData
        
    end% end methods
end % end class def

%% Pst specific thoughts
%{
    Might be a good idea to explore a more dynamic field name approach for 
dealing with the global g structure.
May enable easier inheritance from a base class for get and set data.

Still not really sure about best way to initialize an object.
Requires from old->new and new->old transforms anyway
Kinda seems like 'object' usage just for easier indexing / initializing - and that's fine.

%}

