% Move the file pointer to skip n lines. 
% If current line has not been completely read, it counts as one line
% in the n lines.
% 
% Ningyu Wang, March 2017
function fp = skipLines(fp,n)
for lp1 = 1:1:n
    fgetl(fp) ;
end