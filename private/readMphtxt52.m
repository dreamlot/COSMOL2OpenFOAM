% Read the *.mphtxt file generated by COMSOL. 
% *.mphtxt file contains the mesh inforamtion
% in COMSOL. 
% This program is for single object model. It 
% is intended to replace the stl file. Because
% COMSOL does not export stl files correctly.
% 
% Inputs:
%   file_name, string:      file name, with .mphtxt
% 
% Outputs:
%   V, matrix:              n1-by-3 matrix of nodes 3D coordinates
%   elem_data, struct array:      
%       ID, integer:        ID of current element
%       dimension, integer: dimension of current element
%       type, string:       type name of current element
% 
%       n_nodes, integer:   number of nodes per element
%       n_elem, integer:    number of elements 
%       nodes, matrix:      member node IDs of each element
% 
%       n_entity, integer:	number of geometric entity indices
%       entity, vector:     master geometric entity indix IDs 
%                           of each element
% 
% Ningyu Wang, March 2017
function [V, elem_data] = readMphtxt52( file_name )
%% Open file
fp = fopen(file_name, 'r') ;
%% Nodes
% skip 18 head lines
skipLines(fp,18) ;

% number of nodes
n_V = fscanf(fp, '%d') ; skipLines(fp,1) ;

% node index adjustment
% add this to the node indeces in faces and tets
adj_V = 1 - fscanf(fp,'%d') ; skipLines(fp,1) ;

% skip lines before the data
skipLines(fp,2) ;

% read vetices info
V = fscanf(fp,'%f %f %f', [3, n_V]) ;
V = V' ;

%% Elements
while( true ) 
    tmp_line = fgetl(fp) ;
%     if endsWith(tmp_line,'# number of element types')
    if ~isempty(strfind(tmp_line,'# number of element types'))
        n_types = ...
            sscanf(tmp_line,'%d # number of element types') ;
        break ;
    end
end
elem_data = struct('ID',{},'dimension',{},'type',{} ...
    ,'n_nodes',{},'n_elem',{},'nodes',{} ...
    ,'n_entity',{},'entity',{}) ;

for lp1 = 1:1:n_types
    elem_data(lp1) = readElement52(fp) ;
    elem_data(lp1).nodes = elem_data(lp1).nodes + adj_V ;
    elem_data(lp1).entity = elem_data(lp1).entity + adj_V ;
end


%% Close file
fclose(fp) ;
    