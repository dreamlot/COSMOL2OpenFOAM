% Read one element from mphtxt file
% 
% Input:
%   fp, file handle:        handle to data file
% Output:
%   elem_data, struct array:      
%       ID, integer:            ID of current element
%       dimension, integer:     dimension of current element
%       type, string:           type name of current element
%       
%       n_nodes, integer:       number of nodes per element
%       n_elem, integer:        number of elements 
%       nodes, matrix:          member node IDs of each element
%       
%       n_parameter, integer:   number of parameter values per element
%       n_p_elem, integer:      number of parameters
%       parameters, matrix:
% 
%       n_entity, integer:      number of domains
%       entity, vector:         master geometric entity indix IDs 
%                               of each element
% 
% Ningyu Wang, March 2017
function [elem_data] = readElement41( fp )
%% Read basic information
% Read ID
while( true ) 
    tmp_line = fgetl(fp) ;
%     pause
    if ~isempty(strfind(tmp_line,'# Type #'))
        elem_data.ID = sscanf(tmp_line, '# Type #%d') ;
        disp(tmp_line) ;
        break ;
    end
end
% Read dimension and type name 
while( true ) 
    tmp_line = fgetl(fp) ;
%     if endsWith(tmp_line,'# type name')
    if ~isempty(strfind(tmp_line,'# type name'))
%         [elem_data.dimension, elem_data.type] ...
%             = sscanf(tmp_line,'%d %c3 # type name') ;
%         tmp_line
        tmp = strsplit(tmp_line) ;
        elem_data.dimension = str2double(tmp(1)) ;
        elem_data.type = tmp(2);
        disp(tmp_line) ;
        break ;
    end
end

%% Read nodes and elements
% Read number of nodes per element
while( true ) 
    tmp_line = fgetl(fp) ;
%     if endsWith(tmp_line,'# number of nodes per element')
    if ~isempty(strfind(tmp_line,'# number of nodes per element'))
        [elem_data.n_nodes] ...
            = sscanf(tmp_line,'%d # number of nodes per element') ;
        break ;
    end
end

% Read number of elements
while( true ) 
    tmp_line = fgetl(fp) ;
%     if endsWith(tmp_line,'# number of elements')
    if ~isempty(strfind(tmp_line,'# number of elements'))
        [elem_data.n_elem] ...
            = sscanf(tmp_line,'%d # number of elements') ;
        break ;
    end
end

% Read nodes for each element
while( true ) 
    tmp_line = fgetl(fp) ;
    if ~isempty(strfind(tmp_line,'# Elements'))
        elem_data.nodes = ...
            fscanf(fp,'%d', [elem_data.n_nodes, elem_data.n_elem]) ;
        elem_data.nodes = elem_data.nodes' ;
        disp(tmp_line) ;
        break ;
    end
end

%% Read parameters
% Read number of parameter values per element
while( true ) 
    tmp_line = fgetl(fp) ;
%     if endsWith(tmp_line,'# number of geometric entity indices')
    if ~isempty(strfind(tmp_line,'# number of parameter values per element'))
        [elem_data.n_parameter] ...
            = sscanf(tmp_line,'%d # number of parameter values per element') ;
        disp(tmp_line) ;
        break ;
    end
end

% Read number of parameters
while( true ) 
    tmp_line = fgetl(fp) ;
%     if endsWith(tmp_line,'# number of parameters')
    if ~isempty(strfind(tmp_line,'# number of parameters'))
        [elem_data.n_p_elem] ...
            = sscanf(tmp_line,'%d # number of parameters') ;
        break ;
    end
end

% Read nodes for each element
if elem_data.n_parameter==0 || elem_data.n_p_elem==0
     elem_data.parameters = [] ;
else
    while( true ) 
        tmp_line = fgetl(fp) ;
        if ~isempty(strfind(tmp_line,'# Parameters'))
            elem_data.parameters = ...
                fscanf(fp,'%d', [elem_data.n_parameter, elem_data.n_p_elem]) ;
            elem_data.parameters = elem_data.parameters' ;
            disp(tmp_line) ;
            break ;
        end
    end
end

%% Read domains
% Read number of domains
while( true ) 
    tmp_line = fgetl(fp) ;
%     if endsWith(tmp_line,'# number of domains')
    if ~isempty(strfind(tmp_line,'# number of domains'))
        [elem_data.n_entity] ...
            = sscanf(tmp_line,'%d # number of domains') ;
        break ;
    end
end

% Read nodes for each element
while( true ) 
    tmp_line = fgetl(fp) ;
    if ~isempty(strfind(tmp_line,'# Domains'))
        elem_data.entity = ...
            fscanf(fp,'%d', elem_data.n_entity) ;
        disp(tmp_line) ;
        break ;
    end
end
disp('Finish reading one element type.') ;


    