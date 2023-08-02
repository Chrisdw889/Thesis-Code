%%% Version - v1.00
%%% Author  - Christopher Weaver

%%% Summary:
%%% This class is a Super class for entire class system.

classdef (Abstract) Study
    
    properties (Abstract)
        directory
    end
    
    methods(Abstract)
        
        import(obj)
        
    end
end

