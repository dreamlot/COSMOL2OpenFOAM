% Copy all lines from a file into another file.
% Input & Output:
%   fp_from, file handle:   source file
%   fp_to, file handle:     target file
function [fp_from, fp_to] = copyFileH( fp_from, fp_to )
while( ~feof(fp_from) )
    fprintf(fp_to,fgetl(fp_from)) ;
    fprintf(fp_to,'\n') ;
end