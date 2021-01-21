% Find whether the 4th node is in the positive direction of the face
% defined by the first 3 nodes.
% 
% Input:
%   face_nodes, row vector:     node IDs in target face
%   points, n-by-3 matrix:      node coordinates
%   cell_nodes, row vector:     node IDs in target cell
% output:
%   isPositiveVol, integer:     0: negative, 1: positive
% 
% Ningyu Wang,
% 2017-04
function isPositiveVol = isPositiveVol(face_nodes, points, cell_nodes)
% Get the ID of the 4th node.
lastNode = sum(~ismember(cell_nodes,face_nodes) * cell_nodes') ;

dir_normal = cross( (points(face_nodes(2),:)-points(face_nodes(1),:)) ...
                   ,(points(face_nodes(3),:)-points(face_nodes(2),:)) ) ;
dir_14 = points(lastNode,:) - points(face_nodes(1),:) ;

isPositiveVol = ( dot(dir_normal,dir_14) > 0 ) ;