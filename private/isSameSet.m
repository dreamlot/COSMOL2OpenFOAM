% Check whether two sets are constituted by the same elements.
% 
% Ningyu Wang, March 2017
function isSameSet = isSameSet( set1, set2, n_elem )
switch nargin
    case 2
    isSameSet = min( max(set1==set2(1)) ...
        , min( max(set1==set2(2)), max(set1==set2(3)) ) ) ;
    case 3
        isSameSet = min( max(set1==(set2(1))) ...
            , max(set1==(set2(2))) ) ;
        for lp1 = 2:1:n_elem
            if ~isSameSet
                break ;
            end
            isSameSet = min( isSameSet, max(set1==(set2(lp1))) ) ;
        end
    otherwise
        error('Wrong numbers of input.') ;
end