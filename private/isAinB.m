% Check whether set A is in set B.
% 
% Ningyu
function isAinB = isAinB(A,B)
isAinB = min( ismember(A,B) );