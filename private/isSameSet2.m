% Check whether two sets are constituted by the same elements.
% 
% Ningyu Wang, March 2017
function isSameSet = isSameSet2( set1, set2 )

isSameSet = min(isAinB(set1,set2), isAinB(set2,set1)) ;