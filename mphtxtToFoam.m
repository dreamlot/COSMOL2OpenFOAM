% Convert COMSOL mphtxt to OpenFOAM polyMesh files.
% Tested with COMSOL 4.1, 5.2, 5.3, 5.3a(as 5.3)
% input:
%   path:               Path of 
%   convertToMeters:    Optional. Defalt 1.
%                       Same as in OpenFOAM blockMesh.
%                       If the unit of the file is micron, give 1e-6.
%   COMSOL_version:     Automatically detected. Define to accelerate.
%                       Version of COMSOL file. 4.2 or 5.1
% 
% Ningyu Wang, April 2017
% 
% 2017-04-23-01-10 modify:
%       Use path instead of file_name. The program will automatically find
%       the first file with the extension *.mphtxt.
%       Parallelized.
% 2017-04-10-11:15 modify: 
%       Use face_flag for flagging.
%       Add ordering in'Calculate face-based relations' after 1.
%       Add grouping when ordering. 
%       Modify the 'Set owner and neighbour' to loop in each group.
%       Modify the 'Look up every face' to do similarly as in last step.
% 
% function mphtxtToFoam(file_name, convertToMeters, COMSOL_version)
function mphtxtToFoam(path, convertToMeters, COMSOL_version)
tmp = dir([path,'/*.mphtxt']) ;
file_name = char( [path,'/',tmp(1).name] ) ;     % compatability for Linux
% file_name = '../../run/round_oil_blob3_tst/round_oil_blob3_tst.mphtxt' ;
time0 = cputime() ;
%% Read COMSOL file
if nargin ==1
    convertToMeters = 1 ;
    COMSOL_version = findVersion(file_name) ;
elseif nargin == 2
    COMSOL_version = findVersion(file_name) ;
elseif nargin == 3
    switch COMSOL_version
        case {4.1, 41, 5.2, 52, 5.3, 53} % pre-programmed and tested
        otherwise               % automatic detect file format
            COMSOL_version = findVersion(file_name) ;
    end
else
    error('Too many input arguments.') ;
end

switch COMSOL_version
    case {4.1, 41}
        [points, elem_data] = readMphtxt41(file_name);
    case {5.2, 52}
        [points, elem_data] = readMphtxt52(file_name);
    case {5.3, 53}
        [points, elem_data] = readMphtxt53(file_name);
end

for lp1 = 1:1:numel(elem_data)
    switch char(elem_data(lp1).type)
        case 'vtx'
            ID_vtx = lp1 ;
        case 'edg'
            ID_edg = lp1 ;
        case 'tri'
            ID_tri = lp1 ;
        case 'tet'
            ID_tet = lp1 ;
    end
end
time1 = cputime() ;
disp(['Read mphtxt file "',file_name,'" in ',num2str(time1-time0),'s.']) ;
%% Initialize folder
% Ningyu 2020APR01
oldfolder = cd([path,'/constant']);
system('mkdir polyMesh');
cd(oldfolder);
%% Export file: points
fp = fopen([path,'/constant/polyMesh/points'],'w') ;
copyFile(fp,'private/points') ;
n_points = size(points,1) ;
fprintf(fp,'\n\n') ;
fprintf(fp,'%d\n', n_points) ;
fprintf(fp,'(\n') ;
% format long ;
fprintf(fp,'(%e %e %e)\n',points' * convertToMeters) ;
fprintf(fp,')\n\n\n') ;
fprintf(fp,'// ************************************************************************* //') ;
fclose(fp);
time2 = cputime() ;
disp(['Export "polyMesh/points" in ',num2str(time2-time1),'s.']) ;
time1 = time2 ;
%% Calculate face-based relations
% Get a list of faces in each cell.

% 1. Get a large table of faces per cell.   ========================
% % fprintf(fp,'1. Get a large table of faces per cell.\n') ;
disp('1. Get a large table of faces per cell.');
% cell ID, neighbour, nodes IDs
% neighbour:    < 0 : cell ID - modified : now internal faces
%               = 0 : Not covered. Also boundary. - now boundary
%               > 0 : index of surface * (-1) - not possible now
face_list = zeros(elem_data(ID_tet).n_elem*4,1+1+3) ; 
for lp1 = 1:1:elem_data(ID_tet).n_elem
    face_list(lp1*4-3,:) = [lp1,0,elem_data(4).nodes(lp1,[1,2,3])] ;
    face_list(lp1*4-2,:) = [lp1,0,elem_data(4).nodes(lp1,[2,3,4])] ;
    face_list(lp1*4-1,:) = [lp1,0,elem_data(4).nodes(lp1,[3,4,1])] ;
    face_list(lp1*4  ,:) = [lp1,0,elem_data(4).nodes(lp1,[4,1,2])] ;
end

% 2. Order the nodes for each face. Group the faces by the first node ID.
disp('2. Order and Group') ;
% create groups
face_group = struct('face_list',{}) ; % cell ID, face ID, nodes IDs
for lp1 = 1:1:n_points
    face_group(lp1).face_list = [] ;
end
% ordering and grouping
n_face_list = size(face_list,1) ;
for lp1 = 1:1:n_face_list
    face_list(lp1,3:5) = sort( face_list(lp1,3:5) ) ;
    face_group(face_list(lp1,3)).face_list = ...
        [ face_group(face_list(lp1,3)).face_list; 
          face_list(lp1,1), lp1, face_list(lp1,3:5)                 ] ;
end

% 3. Set owner and neighbour.               ========================
disp('3. Set owner and neighbour.') ;
%   Find boundaries.
%   Delete duplicated faces. 
faces = [] ;
owner = [] ;
neighbour = [] ;
boundary = [] ;

% A tmp flag array. later apply: face_list(face_flag2,2) = 1 ;
face_flag = [];    
% parfor lp1 = 1:1:n_points
for lp1 = 1:1:n_points
    % For each face, compare with all previous faces to check whether it
    % has appeared earlier. If yes, this face is also in another
    % tetrahedron.
%     disp([num2str(lp1),'/',num2str(n_points)]) ;
    
    for lp2 = 1:1:size(face_group(lp1).face_list,1)
        
        for lp3 = 1:1:lp2-1
            % If the 2nd element in a face is not zero, this face is already
            % marked.
            if( face_group(lp1).face_list(lp2,3:5) == face_group(lp1).face_list(lp3,3:5) )
                tmp = [ face_group(lp1).face_list(lp2,2) ;
                        face_group(lp1).face_list(lp3,2) ] ;
                face_flag = [face_flag; tmp] ;

                owner =     [owner;     face_group(lp1).face_list(lp2,1)    ] ;
                neighbour = [neighbour; face_group(lp1).face_list(lp3,1)    ] ;

                % face normal vector points from owner to neighbour                    
                if isPositiveVol(face_group(lp1).face_list(lp2,3:5), points ...
                        , elem_data(ID_tet).nodes(face_group(lp1).face_list(lp2,1),:) )
                    faces =     [faces;     face_group(lp1).face_list(lp2,[3 5 4])  ] ;
                else
                    faces =     [faces;     face_group(lp1).face_list(lp2,3:5)  ] ;
                end
                break ;
            end
        end
    end
end
% face_flag = sort(face_flag) ;
face_list(face_flag,2) = -1 ;

% Now everything left are boundary faces


% 4. Look up every face in the boundary information to know which surface
% it belongs to.                             ========================
disp('4. Look up every face in the boundary information to know which surface it belongs to.') ;
boundary_flag = [] ;
boundary_face_flag = find(face_list(:,2)==0) ;
n_boundary_face_flag = numel(boundary_face_flag) ;

% Order the nodes in each surface face.
for lp1 = 1:1:elem_data(ID_tri).n_elem
    elem_data(ID_tri).nodes(lp1,:) = sort( elem_data(ID_tri).nodes(lp1,:) ) ;
end
parfor lp1_tmp = 1:1:n_boundary_face_flag
    lp1 = boundary_face_flag(lp1_tmp) ;
%     disp([num2str(lp1_tmp),'/',num2str(n_boundary_face_flag)]) ;
    for lp2 = 1:1:elem_data(ID_tri).n_elem
        if ( face_list(lp1,3:5) == elem_data(ID_tri).nodes(lp2,:) )
            tmp = [lp1,elem_data(ID_tri).entity(lp2)];
            boundary_flag = [boundary_flag; tmp] ;
        end
    end
end
face_list(boundary_flag(:,1),2) = boundary_flag(:,2) ;
% 5. Add boundary faces to faces list, owner list, and boundary list.
%                                           ========================
disp('5. Add boundary faces to faces list, owner list, and boundary list.') ;

n_boundary = max(elem_data(ID_tri).entity) ; % number of different sections
nFaces = zeros(n_boundary,1) ;  % used in OpenFOAM format
startFace = nFaces ;            % used in OpenFOAM format
for lp1 = 1:1:n_boundary
    if ~isAinB(lp1,elem_data(ID_tri).entity)
        continue ;
    end
    % get all faces in this boudnary set
    listPerBoundary = find(face_list(:,2)==lp1) ;
    startFace(lp1) = size(faces,1) + 1 ;
    nFaces(lp1) = numel(listPerBoundary) ;
    for lp2 = 1:1:nFaces(lp1)
        owner = [owner; face_list(listPerBoundary(lp2),1)] ;
        if isPositiveVol(face_list(listPerBoundary(lp2),3:5), points ...
                , elem_data(ID_tet).nodes(face_list(listPerBoundary(lp2),1),:))
            faces = [faces; face_list(listPerBoundary(lp2),[3 5 4]) ] ;
        else
            faces = [faces; face_list(listPerBoundary(lp2),3:5) ];
        end
    end
end

% 6. Summarize.                             ========================
disp('6. Summarize.');

n_faces = size(faces,1) ;           % number of all faces
n_owner = n_faces ;                 % number of faces in file 'owner'
n_neighbour = size(neighbour,1) ;   % number of faces in file 'neighbour'
n_internal = n_neighbour ;          % number of internal faces
time2 = cputime() ;
disp(['Calculate face-based relations in ',num2str(time2-time1),'s.']) ;
time1 = time2 ;

%% Export file: faces
fp = fopen([path,'/constant/polyMesh/faces'],'w') ;
copyFile(fp,'private/faces') ;
fprintf(fp,'\n\n') ;
fprintf(fp,'%d\n', n_faces) ;
fprintf(fp,'(\n') ;
fprintf(fp,'3(%d %d %d)\n',faces'-1) ;   % OpenFOAM starts at 0
fprintf(fp,')\n\n\n') ;
fprintf(fp,'// ************************************************************************* //') ;
fclose(fp);
time2 = cputime() ;
disp(['Export "polyMesh/faces" in ',num2str(time2-time1),'s.']) ;
time1 = time2 ;
%% Export file: owner
fp = fopen([path,'/constant/polyMesh/owner'],'w') ;
copyFile(fp,'private/owner') ;
fprintf(fp,'    note        "nPoints: %d nCells: %d nFaces: %d nInternalFaces: %d";\n' ...
    , n_points, elem_data(ID_tet).n_elem, n_faces, n_neighbour) ;
fprintf(fp,'    location    "constant/polyMesh";\n') ;
fprintf(fp,'    object      owner;\n') ;
fprintf(fp,'}\n') ;
fprintf(fp,'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //') ;
fprintf(fp,'\n\n') ;
fprintf(fp,'%d\n', n_owner) ;
fprintf(fp,'(\n') ;
fprintf(fp,'%d\n',owner-1) ;            % OpenFOAM starts at 0
fprintf(fp,')\n\n\n') ;
fprintf(fp,'// ************************************************************************* //') ;
fclose(fp);
time2 = cputime() ;
disp(['Export "polyMesh/owner" in ',num2str(time2-time1),'s.']) ;
time1 = time2 ;
%% Export file: neighbour
fp = fopen([path,'/constant/polyMesh/neighbour'],'w') ;
copyFile(fp,'private/neighbour') ;
fprintf(fp,'    note        "nPoints: %d nCells: %d nFaces: %d nInternalFaces: %d";\n' ...
    , n_points, elem_data(ID_tet).n_elem, n_faces, n_internal) ;
fprintf(fp,'    location    "constant/polyMesh";\n') ;
fprintf(fp,'    object      neighbour;\n') ;
fprintf(fp,'}\n') ;
fprintf(fp,'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //') ;
fprintf(fp,'\n\n') ;
fprintf(fp,'%d\n', n_neighbour) ;
fprintf(fp,'(\n') ;
fprintf(fp,'%d\n',neighbour-1) ;        % OpenFOAM starts at 0
fprintf(fp,')\n\n\n') ;
fprintf(fp,'// ************************************************************************* //') ;
fclose(fp);
time2 = cputime() ;
disp(['Export "polyMesh/neighbour" in ',num2str(time2-time1),'s.']) ;
time1 = time2 ;
%% Export file: boundary
fp = fopen([path,'/constant/polyMesh/boundary'],'w') ;
copyFile(fp,'private/boundary') ;
fprintf(fp,'\n\n') ;
fprintf(fp,'%d\n', n_boundary) ;
fprintf(fp,'(\n') ;
for lp1=1:1:n_boundary
    if ~isAinB(lp1,elem_data(ID_tri).entity)
        continue ;
    end
    fprintf(fp,'    patch%d\n', lp1);
    fprintf(fp,'    {\n') ;
    fprintf(fp,'        type            patch;\n') ;
    fprintf(fp,'        nFaces          %d;\n', nFaces(lp1)) ;
    fprintf(fp,'        startFace       %d;\n', startFace(lp1)-1) ;
    fprintf(fp,'    }\n') ;
end
fprintf(fp,')\n\n\n') ;
fprintf(fp,'// ************************************************************************* //') ;
fclose(fp);
% time2 = cputime() ;
% disp(['Export "polyMesh/boundary" in ',num2str(time2-time1),'s.']) ;
% time1 = time2 ;


% disp(['mphtxtToFoam succeed in ',num2str(time2-time0),'s.']) ;

% Save log file
fp = fopen('log_mphToFoam.txt','a+') ;
% fprintf(fp,'================ %u-%u-%u-%u:%u:%u ================\n', clock) ;
fprintf(fp,'file name = %s\n', file_name) ;
fprintf(fp,'convertToMeters = %f\n', convertToMeters) ;
fprintf(fp,'COMSOL_version = %f\n', COMSOL_version) ;
% fprintf(fp,'Total run time: %f s\n', cputime()-time0 ) ;
fprintf(fp,'nodes: %u, cells: %u, faces: %u, boundary faces: %u\n\n' ...
        , n_points, elem_data(ID_tet).n_elem, n_owner, n_owner-n_neighbour ) ;
fclose(fp) ;
