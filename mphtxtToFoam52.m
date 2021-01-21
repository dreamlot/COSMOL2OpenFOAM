% Convert COMSOL mphtxt to OpenFOAM polyMesh files.
% 
% Ningyu Wang, March 2017
function mphtxtToFoam52(file_name, convertToMeters)
time0 = cputime() ;
%% Read COMSOL file

[points, elem_data] = readMphtxt52(file_name);
for lp1 = 1:1:numel(elem_data)
    switch elem_data(lp1).type
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
system('mkdir polyMesh');
%% Export file: points
fp = fopen('polyMesh/points','w') ;
copyFile(fp,'private/points') ;
n_points = size(points,1) ;
fprintf(fp,'\n\n') ;
fprintf(fp,'%d\n', n_points) ;
fprintf(fp,'(\n') ;
format long ;
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

% cell ID, neighbour, nodes IDs
% neighbour:    > 0 : cell ID
%               = 0 : Not covered. Also boundary.
%               < 0 : index of surface * (-1)
face_list = zeros(elem_data(ID_tet).n_elem*4,1+1+3) ; 
for lp1 = 1:1:elem_data(ID_tet).n_elem
    face_list(lp1*4-3,:) = [lp1,0,elem_data(4).nodes(lp1,[1,2,3])] ;
    face_list(lp1*4-2,:) = [lp1,0,elem_data(4).nodes(lp1,[2,3,4])] ;
    face_list(lp1*4-1,:) = [lp1,0,elem_data(4).nodes(lp1,[3,4,1])] ;
    face_list(lp1*4  ,:) = [lp1,0,elem_data(4).nodes(lp1,[4,1,2])] ;
end

% 2. Set owner and neighbour.               ========================

%   Find boundaries.
%   Delete duplicated faces. 
faces = [] ;
owner = [] ;
neighbour = [] ;
boundary = [] ;
% tmp_boundary = [ elem_data(ID_tri).entity, elem_data(ID_tri).nodes ] ; 
% n_t_b = size(tmp_boundary,1) ;

for lp1 = 1:1:size(face_list,1)
    % For each face, compare with all previous faces to check whether it
    % has appeared earlier. If yes, this face is also in another
    % tetrahedron.
    for lp2 = 1:1:lp1-1
        % If the 2nd element in a face is not zero, this face is already
        % marked.
        if ~face_list(lp2,2)
            if( isAinB( face_list(lp1,3:5), face_list(lp2,3:5) ) )
                face_list(lp2,2) = face_list(lp1,1) ;
                face_list(lp1,2) = 0.1 ; % non-zero, int32(0.1) is zero
                owner =     [owner;     face_list(lp2,1)    ] ;
                neighbour = [neighbour; face_list(lp2,2)    ] ;
                faces =     [faces;     face_list(lp2,3:5)  ] ;
                break ;
            end
        end
    end
end
% Now everything left are boundary faces
tmp_boundary = face_list( (face_list(:,2)==0), : ) ;

% 3. Look up every face in the boundary information to know which surface
%                                           ========================
% it belongs to.
for lp1 = 1:1:size(tmp_boundary,1)
    for lp2 = 1:1:elem_data(ID_tri).n_elem
        if isAinB( tmp_boundary(lp1,3:5), elem_data(ID_tri).nodes(lp2,:) )
            tmp_boundary(lp1,2) = elem_data(ID_tri).entity(lp2) ;
        end
    end
end

% 4. Add boundary faces to faces list, owner list, and boundary list.
%                                           ========================

n_boundary = max(elem_data(ID_tri).entity) ; % number of different sections
nFaces = zeros(n_boundary,1) ;  % used in OpenFOAM format
startFace = nFaces ;            % used in OpenFOAM format
for lp1 = 1:1:n_boundary
    % get all faces in this boudnary set
    tmp = (tmp_boundary(:,2)==lp1) ;
    startFace(lp1) = size(faces,1) + 1 ;
    nFaces(lp1) = sum(tmp) ;
    boundary = [boundary; tmp_boundary(tmp,:)];
    owner = [owner; tmp_boundary(tmp,1)] ;
    faces = [faces; tmp_boundary(tmp,3:5)] ;
end

% 5. Make sure the face ID in 'owner' is lower than in 'neighbour'.
for lp1 = 1:1:numel(n_neighbour)
    if( owner(lp1) > neighbour(lp1) )
        tmp = owner(lp1) ;
        owner(lp1) = neighbour(lp2) ;
        neighbour(lp1) = tmp ;
    end
end

% 6. Summarize.                             ========================

n_faces = size(faces,1) ;           % number of all faces
n_owner = n_faces ;                 % number of faces in file 'owner'
n_neighbour = size(neighbour,1) ;   % number of faces in file 'neighbour'
n_internal = n_neighbour ;          % number of internal faces
time2 = cputime() ;
disp(['Calculate face-based relations in ',num2str(time2-time1),'s.']) ;
time1 = time2 ;
%% Export file: faces
fp = fopen('polyMesh/faces','w') ;
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
fp = fopen('polyMesh/owner','w') ;
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
fp = fopen('polyMesh/neighbour','w') ;
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
fp = fopen('polyMesh/boundary','w') ;
copyFile(fp,'private/boundary') ;
fprintf(fp,'\n\n') ;
fprintf(fp,'%d\n', n_boundary) ;
fprintf(fp,'(\n') ;
for lp1=1:1:n_boundary
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
time2 = cputime() ;
disp(['Export "polyMesh/boundary" in ',num2str(time2-time1),'s.']) ;
time1 = time2 ;