% Copy all lines from a file into another file.
% Input & Output:
%   name_from, string:      source file
%   fp_to, file handle:     target file
function fp_to = copyFile( fp_to, name_from  )
fp_from = fopen(name_from, 'r') ;
while( ~feof(fp_from) )
    fprintf(fp_to,fgetl(fp_from)) ;
    fprintf(fp_to,'\n') ;
end
fclose(fp_from) ;