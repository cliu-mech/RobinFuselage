%% Script that generates meshes of robin fuselage
clc;clear;close all;
%% Function Inputs
meshTypes=["triangular","quadrilateral"]; % mesh type
fileTypes=["plt","vtk"];% output file type (.plt for tecplot, .vtk for paraview)
pylonFlag=true; % flag indicating whether or not pylon exists
coordinateOffset=zeros(3,1); % coordinate of origin in the fuselage frame in global reference frame
L=1.0; % Length of fuselage body
nBodyProfiles=151; % number of longitudinal profiles of fuselage
nAnglesPerProfile=121; % number of angles of discretized points per profile
if(pylonFlag)
    nPylonProfiles=37; % number of longitudinal profiles of pylon
end
%%
pylonLim=[0.4 0.8 1.018];
bodyLim=[0 0.4 0.8 1.9 2.0];
 
% Rows 1, 2, 3 and 4 of the coefficient matrix are for the fuselage
% Rows 5 and 6 are for the pylon
Hcoe=[1.0 -1.0 -0.4 -0.4 1.8 0.0 0.25 1.8;...
    0.0 0.0 0.0 1.0 0.0 0.25 0.0 1.0;...
    1.0 -1.0 -0.8 1.1 1.5 0.05 0.2 0.6;...
    1.0 -1.0 -1.9 0.1 2.0 0.0 0.05 2.0;...
    1.0 -1.0 -0.8 -0.4 3.0 0.0 0.145 3.0;...
    1.0 -1.0 -0.8 0.218 2.0 0.0 0.145 2.0];

Wcoe=[1.0 -1.0 -0.4 -0.4 2.0 0.0 0.25 2.0;...
    0.0 0.0 0.0 1.0 0.0 0.25 0.0 1.0;...
    1.0 -1.0 -0.8 1.1 1.5 0.05 0.2 0.6;...
    1.0 -1.0 -1.9 0.1 2.0 0.0 0.05 2.0;...
    1.0 -1.0 -0.8 -0.4 3.0 0.0 0.166 3.0;...
    1.0 -1.0 -0.8 0.218 2.0 0.0 0.166 2.0];

Zcoe=[1.0 -1.0 -0.4 -0.4 1.8 -0.08 0.08 1.8;...
    0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0;...
    1.0 -1.0 -0.8 1.1 1.5 0.04 -0.04 0.6;...
    0.0 0.0 0.0 1.0 0.0 0.04 0.0 1.0;...
    0.0 0.0 0.0 1.0 0.0 0.125 0.0 1.0;...
    1.0 -1.0 -0.8 1.1 1.5 0.065 0.06 0.6];

Ncoe=[2.0 3.0 0.0 0.4 1.0 0.0 1.0 1.0;...
    0.0 0.0 0.0 1.0 0.0 5.0 0.0 1.0;...
    5.0 -3.0 -0.8 1.1 1.0 0.0 1.0 1.0;...
    0.0 0.0 0.0 1.0 0.0 2.0 0.0 1.0;...
    0.0 0.0 0.0 1.0 0.0 5.0 0.0 1.0;...
    0.0 0.0 0.0 1.0 0.0 5.0 0.0 1.0];

% Fixes from Applied-Scientific-Research/robin-surface-mesh.git
% 1) if there's a 0.0 in the second col, then change the 4th and 5th cols to 1.0
% 2) if there's a 0.0 in C7, change C8 to 1.0, same as above, to prevent nan/inf
% 3) the 0.4..0.8 section (row 2) coefficients in C1 needed to go into C6
% 4) C4 is wrong in the first section of fuse and pyl - it needed to be negative

%% Body
% cosine spacing distribution
angleDist=linspace(180,0,nBodyProfiles);
x_L=cosd(angleDist);
x_L=x_L+ones(size(x_L));
x_L(end)=1.9999;
phi_body=linspace(-360,0,nAnglesPerProfile);

xyz_head_body=zeros(3,1);
xyz_tail_body=zeros(3,1);
x_body=zeros(length(phi_body)-1,length(x_L)-2); % excluding head point and tail point
y_body=x_body;
z_body=x_body;

for idx=1:length(x_L)
    position=find(x_L(idx)>=bodyLim,1,'last');
    % H
    C1=Hcoe(position,1);C2=Hcoe(position,2);
    C3=Hcoe(position,3);C4=Hcoe(position,4);
    C5=Hcoe(position,5);C6=Hcoe(position,6);
    C7=Hcoe(position,7);C8=Hcoe(position,8);
    H=C6+C7*(C1+C2*((x_L(idx)+C3)/C4)^C5)^(1/C8);
    % W
    C1=Wcoe(position,1);C2=Wcoe(position,2);
    C3=Wcoe(position,3);C4=Wcoe(position,4);
    C5=Wcoe(position,5);C6=Wcoe(position,6);
    C7=Wcoe(position,7);C8=Wcoe(position,8);
    W=C6+C7*(C1+C2*((x_L(idx)+C3)/C4)^C5)^(1/C8);
    % Z
    C1=Zcoe(position,1);C2=Zcoe(position,2);
    C3=Zcoe(position,3);C4=Zcoe(position,4);
    C5=Zcoe(position,5);C6=Zcoe(position,6);
    C7=Zcoe(position,7);C8=Zcoe(position,8);
    Z=C6+C7*(C1+C2*((x_L(idx)+C3)/C4)^C5)^(1/C8);
    % N
    C1=Ncoe(position,1);C2=Ncoe(position,2);
    C3=Ncoe(position,3);C4=Ncoe(position,4);
    C5=Ncoe(position,5);C6=Ncoe(position,6);
    C7=Ncoe(position,7);C8=Ncoe(position,8);
    N=C6+C7*(C1+C2*((x_L(idx)+C3)/C4)^C5)^(1/C8);
    
    if(idx==1)% Head point of fuselage
        x_L(idx)=0.0;
        xyz_head_body=[0.0;0.0;Z*L];
        continue;
    elseif(idx==length(x_L))% Tail point of fuselage
        x_L(idx)=2.0;
        xyz_tail_body=[x_L(idx)*L;0.0;Z*L];
        break;
    else% Body points of fuselage
        for iPhi=1:length(phi_body)-1          
            r=((0.25*H*W)^N)^(1/N)/(abs(0.5*H*sind(phi_body(iPhi)))^N+abs(0.5*W*cosd(phi_body(iPhi)))^N)^(1/N);
            y_L=r*sind(phi_body(iPhi));
            z_L=r*cosd(phi_body(iPhi))+Z;                        
            x_body(iPhi,idx-1)=x_L(idx)*L;
            y_body(iPhi,idx-1)=y_L*L;
            z_body(iPhi,idx-1)=z_L*L;
        end
    end
end

%% Pylon
if(pylonFlag)
    dx=1e-6;
    phi_pylon=linspace(-90,90,ceil(nAnglesPerProfile/2));% only half revolution needed (on the top of fuselage)
    %phi_pylon=linspace(-180,180,nLongitudinalAngles);% intact revolution (bottom surface overlaps on fuselage)
    % opt 1: longitudinal profiles that align with fuselage profiles
    %pylonIndices=(x_L>pylonLim(1)&x_L<pylonLim(end));
    %x_L_pylon=[pylonLim(1)+dx x_L(pylonIndices) pylonLim(end)-dx];
    % opt 2: misalignment between pylon profiles and fuselage profiles
    diffProfileLim=pylonLim(end)-pylonLim(1);
    angleDist=linspace(180,0,nPylonProfiles);
    x_L_pylon=cosd(angleDist)*diffProfileLim/2;
    x_L_pylon=x_L_pylon+ones(size(x_L_pylon))*(diffProfileLim/2+pylonLim(1));
    x_L_pylon(end)=x_L_pylon(end)-dx;
    %
    
    nPylonStations=length(x_L_pylon); % caluclate number of stations within pylon limit
    xyz_head_pylon=zeros(3,1);
    xyz_tail_pylon=zeros(3,1);
    x_pylon=zeros(length(phi_pylon),nPylonStations);
    y_pylon=x_pylon;
    z_pylon=x_pylon;

    for idx=1:length(x_L_pylon)
        position=find(x_L_pylon(idx)>=pylonLim, 1, 'last')+4;
        % H
        C1=Hcoe(position,1);C2=Hcoe(position,2);
        C3=Hcoe(position,3);C4=Hcoe(position,4);
        C5=Hcoe(position,5);C6=Hcoe(position,6);
        C7=Hcoe(position,7);C8=Hcoe(position,8);
        H=C6+C7*(C1+C2*((x_L_pylon(idx)+C3)/C4)^C5)^(1/C8);
        % W
        C1=Wcoe(position,1);C2=Wcoe(position,2);
        C3=Wcoe(position,3);C4=Wcoe(position,4);
        C5=Wcoe(position,5);C6=Wcoe(position,6);
        C7=Wcoe(position,7);C8=Wcoe(position,8);
        W=C6+C7*(C1+C2*((x_L_pylon(idx)+C3)/C4)^C5)^(1/C8);
        % Z
        C1=Zcoe(position,1);C2=Zcoe(position,2);
        C3=Zcoe(position,3);C4=Zcoe(position,4);
        C5=Zcoe(position,5);C6=Zcoe(position,6);
        C7=Zcoe(position,7);C8=Zcoe(position,8);
        Z=C6+C7*(C1+C2*((x_L_pylon(idx)+C3)/C4)^C5)^(1/C8);
        % N
        C1=Ncoe(position,1);C2=Ncoe(position,2);
        C3=Ncoe(position,3);C4=Ncoe(position,4);
        C5=Ncoe(position,5);C6=Ncoe(position,6);
        C7=Ncoe(position,7);C8=Ncoe(position,8);
        N=C6+C7*(C1+C2*((x_L_pylon(idx)+C3)/C4)^C5)^(1/C8);
    
        if(idx==1)            
            %xyz_head_pylon=[x_L_pylon(idx);0.0;Z*L];
            xyz_head_pylon=[pylonLim(1)*L;0.0;Z*L];
        end

        if(idx==length(x_L_pylon))
            %xyz_tail_pylon=[x_L_pylon(idx)*L;0.0;Z*L];
            xyz_tail_pylon=[pylonLim(end)*L;0.0;Z*L];
        end

        for iPhi=1:length(phi_pylon)           
            r=((0.25*H*W)^N)^(1/N)/(abs(0.5*H*sind(phi_pylon(iPhi)))^N+abs(0.5*W*cosd(phi_pylon(iPhi)))^N)^(1/N);
            y_L=r*sind(phi_pylon(iPhi));
            z_L=r*cosd(phi_pylon(iPhi))+Z;       
            x_pylon(iPhi,idx)=x_L_pylon(idx)*L;
            y_pylon(iPhi,idx)=y_L*L;
            z_pylon(iPhi,idx)=z_L*L;
        end
    end
end

%% Write meshes into file
WriteVTK_NoPylon(xyz_head_body,xyz_tail_body,x_body,y_body,z_body,2);
WritePLT_NoPylon(xyz_head_body,xyz_tail_body,x_body,y_body,z_body,2);

%% Local Functions
% Generate body meshes
function [triangularBodyPanels,quadrilateralBodyPanels]=GenerateBodyMeshesNoPylon(bodyHead,bodyTail,x_body,y_body,z_body)

nProfiles=size(x_body,2);
nPtsPerProfile=size(x_body,1);
nQuadrilaterPanels=(nProfiles-1)*(nPtsPerProfile-1);
nTriangularPanels=2*nPtsPerProfile;
triangularBodyPanels=zeros(3,3,nTriangularPanels);
quadrilateralBodyPanels=zeros(3,4,nQuadrilaterPanels);

% front
iPanel=1;
for iPt=1:nPtsPerProfile
    nextIdx=mod(iPt,nPtsPerProfile)+1;
    triangularBodyPanels(:,1,iPanel)=bodyHead;
    triangularBodyPanels(:,2,iPanel)=[x_body(iPt,1);y_body(iPt,1);z_body(iPt,1)];    
    triangularBodyPanels(:,3,iPanel)=[x_body(nextIdx,1);y_body(nextIdx,1);z_body(nextIdx,1)];
    iPanel=iPanel+1;
end
% rear
for iPt=1:nPtsPerProfile
    nextIdx=mod(iPt,nPtsPerProfile)+1;
    triangularBodyPanels(:,1,iPanel)=bodyTail;
    triangularBodyPanels(:,2,iPanel)=[x_body(iPt,end);y_body(iPt,end);z_body(iPt,end)];    
    triangularBodyPanels(:,3,iPanel)=[x_body(nextIdx,end);y_body(nextIdx,end);z_body(nextIdx,end)];
    iPanel=iPanel+1;
end
% body
iPanel=1;
for iProfile=1:nProfiles-1
    for iPt=1:nPtsPerProfile
        nextIdx=mod(iPt,nPtsPerProfile)+1;
        quadrilateralBodyPanels(:,1,iPanel)=[x_body(iPt,iProfile);y_body(iPt,iProfile);z_body(iPt,iProfile)];        
        quadrilateralBodyPanels(:,2,iPanel)=[x_body(nextIdx,iProfile);y_body(nextIdx,iProfile);z_body(nextIdx,iProfile)];
        quadrilateralBodyPanels(:,3,iPanel)=[x_body(nextIdx,iProfile+1);y_body(nextIdx,iProfile+1);z_body(nextIdx,iProfile+1)];
        quadrilateralBodyPanels(:,4,iPanel)=[x_body(iPt,iProfile+1);y_body(iPt,iProfile+1);z_body(iPt,iProfile+1)];
        iPanel=iPanel+1;
    end
end
end

% Generate pylon meshes (intact circle)
function [triangularPylonPanels,quadrilateralPylonPanels]=GeneratePylonMeshes(pylonHead,pylonTail,x_pylon,y_pylon,z_pylon)

nProfiles=size(x_pylon,2);
nPtsPerProfile=size(x_pylon,1);
%nPointsAll=2+nProfiles*nPtsPerProfile;% 2 points refer to pylon head and tail points
nSurfaces=(nProfiles-1)*(nPtsPerProfile-1);% surfaces without front and rear parts

nTriangularPanels=2*(nPtsPerProfile-1);
triangularPylonPanels=zeros(3,3,nTriangularPanels);
quadrilateralPylonPanels=zeros(3,4,nSurfaces);

% front
iPanel=1;
for iPoint=1:nPtsPerProfile
    nextIdx=mod(iPoint,nPtsPerProfile)+1;
    triangularPylonPanels(:,1,iPanel)=pylonHead;
    triangularPylonPanels(:,2,iPanel)=[x_pylon(iPoint,1);y_pylon(iPoint,1);z_pylon(iPoint,1)];
    triangularPylonPanels(:,3,iPanel)=[x_pylon(nextIdx,1);y_pylon(nextIdx,1);z_pylon(nextIdx,1)];
    iPanel=iPanel+1;
end
% rear
for iPoint=1:nPtsPerProfile
    nextIdx=mod(iPoint,nPtsPerProfile)+1;
    triangularPylonPanels(:,1,iPanel)=pylonTail;
    triangularPylonPanels(:,2,iPanel)=[x_pylon(iPoint,end);y_pylon(iPoint,end);z_pylon(iPoint,end)];
    triangularPylonPanels(:,3,iPanel)=[x_pylon(nextIdx,end);y_pylon(nextIdx,end);z_pylon(nextIdx,end)];
    iPanel=iPanel+1;
end
% pylon body
iPanel=1;
for iProfile=1:nProfiles-1
    for iPt=1:nPtsPerProfile
        nextIdx=mod(iPt,nPtsPerProfile)+1;
        quadrilateralPylonPanels(:,1,iPanel)=[x_pylon(iPt,iProfile);y_pylon(iPt,iProfile);z_pylon(iPt,iProfile)];
        quadrilateralPylonPanels(:,2,iPanel)=[x_pylon(nextIdx,iProfile);y_pylon(nextIdx,iProfile);z_pylon(nextIdx,iProfile)];
        quadrilateralPylonPanels(:,3,iPanel)=[x_pylon(nextIdx,iProfile+1);y_pylon(nextIdx,iProfile+1);z_pylon(nextIdx,iProfile+1)];
        quadrilateralPylonPanels(:,4,iPanel)=[x_pylon(iPt,iProfile+1);y_pylon(iPt,iProfile+1);z_pylon(iPt,iProfile+1)];
        iPanel=iPanel+1;
    end
end
end

% Generate pylon meshes (half circle)
function [triangularPylonPanels,quadrilateralPylonPanels]=GeneratePylonMeshesHalf(pylonHead,pylonTail,x_pylon,y_pylon,z_pylon)

nProfiles=size(x_pylon,2);
nPtsPerProfile=size(x_pylon,1);
%nPointsAll=2+nProfiles*nPtsPerProfile;% 2 points refer to pylon head and tail points
nSurfaces=(nProfiles-1)*(nPtsPerProfile-1);% surfaces without front and rear parts

nTriangularPanels=2*(nPtsPerProfile-1);
triangularPylonPanels=zeros(3,3,nTriangularPanels);
quadrilateralPylonPanels=zeros(3,4,nSurfaces);

% front
iPanel=1;
for iPoint=1:nPtsPerProfile-1
    nextIdx=iPoint+1;
    triangularPylonPanels(:,1,iPanel)=pylonHead;
    triangularPylonPanels(:,2,iPanel)=[x_pylon(iPoint,1);y_pylon(iPoint,1);z_pylon(iPoint,1)];
    triangularPylonPanels(:,3,iPanel)=[x_pylon(nextIdx,1);y_pylon(nextIdx,1);z_pylon(nextIdx,1)];
    iPanel=iPanel+1;
end
% rear
for iPoint=1:nPtsPerProfile-1
    nextIdx=iPoint+1;
    triangularPylonPanels(:,1,iPanel)=pylonTail;
    triangularPylonPanels(:,2,iPanel)=[x_pylon(iPoint,end);y_pylon(iPoint,end);z_pylon(iPoint,end)];
    triangularPylonPanels(:,3,iPanel)=[x_pylon(nextIdx,end);y_pylon(nextIdx,end);z_pylon(nextIdx,end)];
    iPanel=iPanel+1;
end
% pylon body
iPanel=1;
for iProfile=1:nProfiles-1
    for iPoint=1:nPtsPerProfile-1
        nextIdx=iPoint+1;
        quadrilateralPylonPanels(:,1,iPanel)=[x_pylon(iPoint,iProfile);y_pylon(iPoint,iProfile);z_pylon(iPoint,iProfile)];
        quadrilateralPylonPanels(:,2,iPanel)=[x_pylon(nextIdx,iProfile);y_pylon(nextIdx,iProfile);z_pylon(nextIdx,iProfile)];
        quadrilateralPylonPanels(:,3,iPanel)=[x_pylon(nextIdx,iProfile+1);y_pylon(nextIdx,iProfile+1);z_pylon(nextIdx,iProfile+1)];
        quadrilateralPylonPanels(:,4,iPanel)=[x_pylon(iPoint,iProfile+1);y_pylon(iPoint,iProfile+1);z_pylon(iPoint,iProfile+1)];
        iPanel=iPanel+1;
    end
end
end

%% Functions that write files
function WriteVTK_NoPylon(bodyHead,bodyTail,x_body,y_body,z_body,meshType)
% Specify file name
fileName='RobinFuselage';
% Create and open file
fileID=fopen(strcat(fileName,'.vtk'),'w');
% Part 1: Header
fprintf(fileID,'# vtk DataFile Version 2.0\n');
% Part 2: Title (256 characters maximum, terminated with newline \n character)
fprintf(fileID,'Meshes of robin fuselage.\n');
% Part 3: Data type, either ASCII or BINARY
fprintf(fileID,'ASCII\n');
% Part 4: Geometry/topology. Type is one of: STRUCTURED_POINTS,
% STRUCTURED_GRID, UNSTRUCTURED_GRID, POLYDATA, RECTILINEAR_GRID, AND FIELD
fprintf(fileID,'DATASET POLYDATA\n');
% Part 5: Dataset attributes
%% Write points
nProfiles=size(x_body,2);
nPointsPerProfile=size(x_body,1);
nPointsAll=nProfiles*nPointsPerProfile+2;
fprintf(fileID,'POINTS %d float\n',nPointsAll);
fprintf(fileID,'%f %f %f\n',bodyHead(1),bodyHead(2),bodyHead(3)); % Body head
% Body points
for iProfile=1:nProfiles
    for iPoint=1:nPointsPerProfile
        fprintf(fileID,'%f %f %f\n',x_body(iPoint,iProfile),y_body(iPoint,iProfile),z_body(iPoint,iProfile));
    end
end
fprintf(fileID,'%f %f %f\n',bodyTail(1),bodyTail(2),bodyTail(3)); % Body tail
fprintf(fileID,'\n');

%% Write polygons
nTriangles=nPointsPerProfile*2;
nQuadrilaterals=(nProfiles-1)*nPointsPerProfile;

% Triangles (front)
nTriangleSides=3;
nQuadSides=4;

switch(meshType)
    case 1
        nPolygons=nTriangles+nQuadrilaterals*2;
        Size=(nTriangleSides+1)*nPolygons;
    case 2
        nPolygons=nTriangles+nQuadrilaterals;
        Size=(nTriangleSides+1)*nTriangles+(nQuadSides+1)*nQuadrilaterals;
end

fprintf(fileID,'POLYGONS %d %d\n',nPolygons,Size);
for iPoint=1:nPointsPerProfile
    nextIdx=mod(iPoint,nPointsPerProfile)+1;
    fprintf(fileID,'%d %d %d %d\n',nTriangleSides,0,iPoint,nextIdx);
end

%  Triangles (rear)
idxOffset=nPointsPerProfile*(nProfiles-1)+1;
for iPoint=1:nPointsPerProfile
    nextIdx=mod(iPoint,nPointsPerProfile)+1;
    fprintf(fileID,'%d %d %d %d\n',nTriangleSides,nPointsAll-1,iPoint+idxOffset-1,nextIdx+idxOffset-1);
end

switch(meshType)
    case 1
        % Triangles
        for iProfile=1:nProfiles-1
            idxOffset=(iProfile-1)*nPointsPerProfile;
            for iPoint=1:nPointsPerProfile
                nextIdx=mod(iPoint,nPointsPerProfile)+1;
                fprintf(fileID,'%d %d %d %d\n',nTriangleSides,iPoint+idxOffset,nextIdx+idxOffset,nextIdx+idxOffset+nPointsPerProfile);
                fprintf(fileID,'%d %d %d %d\n',nTriangleSides,iPoint+idxOffset,iPoint+idxOffset+nPointsPerProfile,nextIdx+idxOffset+nPointsPerProfile);
            end
        end
    case 2
        % Quadrilaterals
        for iProfile=1:nProfiles-1
            idxOffset=(iProfile-1)*nPointsPerProfile;
            for iPoint=1:nPointsPerProfile
                nextIdx=mod(iPoint,nPointsPerProfile)+1;
                fprintf(fileID,'%d %d %d %d %d\n',nQuadSides,iPoint+idxOffset,nextIdx+idxOffset,...
                    nextIdx+idxOffset+nPointsPerProfile,iPoint+idxOffset+nPointsPerProfile);
            end
        end    
end

% Close file
fclose(fileID);
end

function WritePLT_NoPylon(bodyHead,bodyTail,x_body,y_body,z_body,meshType)

nProfiles=size(x_body,2);
nPointsPerProfile=size(x_body,1);

% Specify file name
fileName='RobinFuselage';
% Create and open file
fileID=fopen(strcat(fileName,'.plt'),'w');
fprintf(fileID,'TITLE = "ROBIN Fuselage"\n');
fprintf(fileID,'VARIABLES = "x" "y" "z"\n');

%% WRITE ELEMENTS
%% head points (triangles)
fprintf(fileID,'ZONE T="FUSELAGEHEAD" N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n',nPointsPerProfile+1,nPointsPerProfile);
fprintf(fileID,'%f %f %f\n',bodyHead(1),bodyHead(2),bodyHead(3));
for iPoint=1:nPointsPerProfile
    fprintf(fileID,'%f %f %f\n',x_body(iPoint,1),y_body(iPoint,1),z_body(iPoint,1));
end
% Indices of nodes of head panels
for iPoint=1:nPointsPerProfile
    nextIdx=mod(iPoint,nPointsPerProfile)+1;
    fprintf(fileID,'%d %d %d\n',1,iPoint+1,nextIdx+1);
end

%% body points (quadrilaterals)
nBodyPoints=nProfiles*nPointsPerProfile;
switch(meshType)
    case 1
        % Triangles
        nBodyPanels=(nProfiles-1)*nPointsPerProfile*2;
        fprintf(fileID,'ZONE T="FUSELAGEBODY" N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n',nBodyPoints,nBodyPanels);
        for iProfile=1:nProfiles    
            for iPoint=1:nPointsPerProfile
                fprintf(fileID,'%f %f %f\n',x_body(iPoint,iProfile),y_body(iPoint,iProfile),z_body(iPoint,iProfile));
            end
        end
        % Indices of nodes of body panels
        for iProfile=1:nProfiles-1
            idxOffset=(iProfile-1)*nPointsPerProfile;
            for iPoint=1:nPointsPerProfile
                nextIdx=mod(iPoint,nPointsPerProfile)+1;
                fprintf(fileID,'%d %d %d\n',iPoint+idxOffset,nextIdx+idxOffset,nextIdx+idxOffset+nPointsPerProfile);
                fprintf(fileID,'%d %d %d\n',iPoint+idxOffset,nextIdx+idxOffset+nPointsPerProfile,iPoint+idxOffset+nPointsPerProfile);
            end
        end
    case 2
        % Quadrilaterals
        nBodyPanels=(nProfiles-1)*nPointsPerProfile;
        fprintf(fileID,'ZONE T="FUSELAGEBODY" N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n',nBodyPoints,nBodyPanels);
        for iProfile=1:nProfiles    
            for iPoint=1:nPointsPerProfile
                fprintf(fileID,'%f %f %f\n',x_body(iPoint,iProfile),y_body(iPoint,iProfile),z_body(iPoint,iProfile));
            end
        end
        % Indices of nodes of body panels
        for iProfile=1:nProfiles-1
            idxOffset=(iProfile-1)*nPointsPerProfile;
            for iPoint=1:nPointsPerProfile
                nextIdx=mod(iPoint,nPointsPerProfile)+1;
                fprintf(fileID,'%d %d %d %d\n',iPoint+idxOffset,nextIdx+idxOffset,...
                    nextIdx+idxOffset+nPointsPerProfile,iPoint+idxOffset+nPointsPerProfile);
            end
        end
end

%% tail points (triangles)
fprintf(fileID,'ZONE T="FUSELAGETAIL" N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n',nPointsPerProfile+1,nPointsPerProfile);
for iPoint=1:nPointsPerProfile
    fprintf(fileID,'%f %f %f\n',x_body(iPoint,nProfiles),y_body(iPoint,nProfiles),z_body(iPoint,nProfiles));
end
fprintf(fileID,'%f %f %f\n',bodyTail(1),bodyTail(2),bodyTail(3));
% Indices of nodes of tail panels
for iPoint=1:nPointsPerProfile
    nextIdx=mod(iPoint,nPointsPerProfile)+1;
    fprintf(fileID,'%d %d %d\n',nPointsPerProfile+1,iPoint,nextIdx);
end

% Close file
fclose(fileID);

end