%% Script that generates meshes of robin fuselage
clc;clear;close all;
%% Function Inputs
meshTypes=["triangular","quadrilateral"]; % mesh type
fileTypes=["plt","vtk"];% output file type (.plt for tecplot, .vtk for paraview)
pylonFlag=true; % flag indicating whether or not pylon exists
coordinateOffset=zeros(3,1); % coordinate of origin in the fuselage frame in global reference frame
L=1.0; % Length of fuselage body
nPylonProfiles=37; % number of longitudinal profiles of pylon
nBody1stProfiles=37;
nBody2ndProfiles=nPylonProfiles;
nBody3rdProfiles=61;
nBodyProfiles=nBody1stProfiles+nBody2ndProfiles+nBody3rdProfiles-2; % number of longitudinal profiles of fuselage

nPointsPerPylonProfile=37; % number of angles of discretized points per profile
nPointsPerBodyProfile=121;

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

%% Pylon
dx=1e-6;
phi_pylon=linspace(-90,90,nPointsPerPylonProfile);% only half revolution needed (on the top of fuselage)
diffProfileLim=pylonLim(end)-pylonLim(1);
angleDist=linspace(180,0,nPylonProfiles);
x_L_pylon=cosd(angleDist)*diffProfileLim/2;
x_L_pylon=x_L_pylon+ones(size(x_L_pylon))*(diffProfileLim/2+pylonLim(1));
x_L_pylon(end)=x_L_pylon(end)-dx;
%
nPylonStations=length(x_L_pylon); % caluclate number of stations within pylon limit
pylonHead=zeros(3,1);
pylonTail=zeros(3,1);
x_pylon=zeros(length(phi_pylon),nPylonStations-2);
y_pylon=x_pylon;
z_pylon=x_pylon;

for idx=1:length(x_L_pylon)
    position=find(x_L_pylon(idx)>=pylonLim, 1, 'last')+4; % pylon index in matrix starts with the 5th row
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
        pylonHead=[pylonLim(1)*L;0.0;Z*L];
        continue;
    end

    if(idx==length(x_L_pylon))
        %xyz_tail_pylon=[x_L_pylon(idx)*L;0.0;Z*L];
        pylonTail=[pylonLim(end)*L;0.0;Z*L];
        break;
    end

    for iPhi=1:length(phi_pylon)           
        r=((0.25*H*W)^N)^(1/N)/(abs(0.5*H*sind(phi_pylon(iPhi)))^N+abs(0.5*W*cosd(phi_pylon(iPhi)))^N)^(1/N);
        y_L=r*sind(phi_pylon(iPhi));
        z_L=r*cosd(phi_pylon(iPhi))+Z;       

        x_pylon(iPhi,idx-1)=x_L_pylon(idx)*L;
        y_pylon(iPhi,idx-1)=y_L*L;
        z_pylon(iPhi,idx-1)=z_L*L;
    end
end

%% Body
x_L_body1=linspace(180,90,nBody1stProfiles);
x_L_body2=x_L_pylon;
x_L_body2(end)=pylonLim(end);
x_L_body3=linspace(90,0,nBody3rdProfiles);

x_L_body1=cosd(x_L_body1)*pylonLim(1);
x_L_body1=x_L_body1+ones(size(x_L_body1))*pylonLim(1);
diffLim=bodyLim(end)-pylonLim(end);
x_L_body3=cosd(x_L_body3)*diffLim;
x_L_body3=x_L_body3+ones(size(x_L_body3))*pylonLim(end);
x_L_body3(end)=1.9999;
x_L=[x_L_body1 x_L_body2(2:end-1) x_L_body3];

bodyHead=zeros(3,1);
bodyTail=zeros(3,1);

x_body=zeros(nPointsPerBodyProfile,nBodyProfiles-2);
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

    if(idx==1)
        bodyHead=[0.0;0.0;Z*L];
        continue;
    elseif(idx==length(x_L))
        x_L(idx)=2.0;
        bodyTail=[x_L(idx)*L;0.0;Z*L];
        break;
    else
        if(x_L(idx)<=pylonLim(1)||x_L(idx)>=pylonLim(end))
            phi_body=linspace(0,360,nPointsPerBodyProfile+1);
            for iPhi=1:length(phi_body)        
                r=((0.25*H*W)^N)^(1/N)/(abs(0.5*H*sind(phi_body(iPhi)))^N+abs(0.5*W*cosd(phi_body(iPhi)))^N)^(1/N);
                y_L=r*sind(phi_body(iPhi));
                z_L=r*cosd(phi_body(iPhi))+Z;                        
                x_body(iPhi,idx-1)=x_L(idx)*L;
                y_body(iPhi,idx-1)=y_L*L;
                z_body(iPhi,idx-1)=z_L*L;
            end
        else
            idxPylon=idx-length(x_L_body1);
            yBound=y_pylon(end,idxPylon)/L;
            zBound=z_pylon(end,idxPylon)/L-Z;
            startAngle=atan2d(yBound,zBound);
            endAngle=360-startAngle;
            phi_body=linspace(startAngle,endAngle,nPointsPerBodyProfile+1);

            x_body(1,idx-1)=x_pylon(end,idxPylon);
            y_body(1,idx-1)=y_pylon(end,idxPylon);
            z_body(1,idx-1)=z_pylon(end,idxPylon);
            x_body(end,idx-1)=x_pylon(1,idxPylon);
            y_body(end,idx-1)=y_pylon(1,idxPylon);
            z_body(end,idx-1)=z_pylon(1,idxPylon);
            for iPhi=2:length(phi_body)-1         
                r=((0.25*H*W)^N)^(1/N)/(abs(0.5*H*sind(phi_body(iPhi)))^N+abs(0.5*W*cosd(phi_body(iPhi)))^N)^(1/N);
                y_L=r*sind(phi_body(iPhi));
                z_L=r*cosd(phi_body(iPhi))+Z;                        
                x_body(iPhi,idx-1)=x_L(idx)*L;
                y_body(iPhi,idx-1)=y_L*L;
                z_body(iPhi,idx-1)=z_L*L;
            end
        end        
    end
end

%%
WriteVTK(2,bodyHead,bodyTail,x_body,y_body,z_body,pylonHead,pylonTail,x_pylon,y_pylon,z_pylon);
WritePLT(2,bodyHead,bodyTail,x_body,y_body,z_body,pylonHead,pylonTail,x_pylon,y_pylon,z_pylon);

%% Local functions that generate .VTK and .PLT files
%% Write VTK File
function WriteVTK(meshType,bodyHead,bodyTail,x_body,y_body,z_body,pylonHead,pylonTail,x_pylon,y_pylon,z_pylon)
idxPoint=0;
% Specify file name
fileName='RobinFuselagePylon';
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
% Calculate number of points and elements
nBodyProfiles=size(x_body,2);
nPointsPerBodyProfile=size(x_body,1);
nBodyPointsAll=2+nBodyProfiles*nPointsPerBodyProfile;
nPylonProfiles=size(x_pylon,2);
nPointsPerPylonProfile=size(x_pylon,1);
nPylonPointsAll=2+nPylonProfiles*nPointsPerPylonProfile;
nPointsAll=nBodyPointsAll+nPylonPointsAll;
% Side numbers of polygons
nPylonTriangles=2*(nPointsPerPylonProfile-1);
nBodyTriangles=2*(nPointsPerBodyProfile-1);
nTriangles=nPylonTriangles+nBodyTriangles;
nPylonQuadrilaterals=(nPointsPerPylonProfile-1)*(nPylonProfiles-1);
nBodyQuadrilaterals=(nPointsPerBodyProfile-1)*(nBodyProfiles-1);
nQuadrilaterals=nPylonQuadrilaterals+nBodyQuadrilaterals;
nTriangleSides=3;
nQuadSides=4;
switch(meshType)
    case 1
        % triangular elements
        nPolygons=nTriangles+nQuadrilaterals*2;
        Size=(nTriangleSides+1)*nPolygons;
    case 2
        % quadrilateral elements
        nPolygons=nTriangles+nQuadrilaterals;
        Size=(nTriangleSides+1)*nTriangles+(nQuadSides+1)*nQuadrilaterals;
end
fprintf(fileID,'POINTS %d float\n',nPointsAll);
%%%%%%%%%%%%%%%%% Body points %%%%%%%%%%%%%%%%%%%%
fprintf(fileID,'%f %f %f\n',bodyHead(1),bodyHead(2),bodyHead(3)); % Body head
idxBodyHead=idxPoint; % Body head index
idxPoint=idxPoint+1;
for iProfile=1:nBodyProfiles
    for iPoint=1:nPointsPerBodyProfile
        fprintf(fileID,'%f %f %f\n',x_body(iPoint,iProfile),y_body(iPoint,iProfile),z_body(iPoint,iProfile));
        idxPoint=idxPoint+1;
    end
end
fprintf(fileID,'%f %f %f\n',bodyTail(1),bodyTail(2),bodyTail(3)); % Body tail
idxBodyTail=idxPoint; % Body tail index
idxPoint=idxPoint+1;
%%%%%%%%%%%%%%%%%% Pylon points %%%%%%%%%%%%%%%%%%%
fprintf(fileID,'%f %f %f\n',pylonHead(1),pylonHead(2),pylonHead(3)); % Pylon head
idxPylonHead=idxPoint;
idxPoint=idxPoint+1;
for iProfile=1:nPylonProfiles
    for iPoint=1:nPointsPerPylonProfile
        fprintf(fileID,'%f %f %f\n',x_pylon(iPoint,iProfile),y_pylon(iPoint,iProfile),z_pylon(iPoint,iProfile));
        idxPoint=idxPoint+1;
    end
end
fprintf(fileID,'%f %f %f\n',pylonTail(1),pylonTail(2),pylonTail(3)); % Pylon tail
idxPylonTail=idxPoint;
fprintf(fileID,'\n');
%% Write polygons
fprintf(fileID,'POLYGONS %d %d\n',nPolygons,Size);
%%%%%%%%%%%%%%%%%%% Body polygons %%%%%%%%%%%%%%%%%%%%%%%%
% Triangles (front)
for iPoint=1:nPointsPerBodyProfile-1
    nextIdx=iPoint+1;
    fprintf(fileID,'%d %d %d %d\n',nTriangleSides,idxBodyHead,iPoint,nextIdx);
end
%  Triangles (rear)
idxOffset=nPointsPerBodyProfile*(nBodyProfiles-1);
for iPoint=1:nPointsPerBodyProfile-1
    nextIdx=iPoint+1;
    fprintf(fileID,'%d %d %d %d\n',nTriangleSides,idxBodyTail,iPoint+idxOffset,nextIdx+idxOffset);
end
% Body polygons (depends on mesh type)
switch(meshType)
    case 1
        % Triangles
        for iProfile=1:nBodyProfiles-1
            idxOffset=(iProfile-1)*nPointsPerBodyProfile;
            for iPoint=1:nPointsPerBodyProfile-1
                nextIdx=iPoint+1;
                fprintf(fileID,'%d %d %d %d\n',nTriangleSides,iPoint+idxOffset,nextIdx+idxOffset,nextIdx+idxOffset+nPointsPerBodyProfile);
                fprintf(fileID,'%d %d %d %d\n',nTriangleSides,iPoint+idxOffset,iPoint+idxOffset+nPointsPerBodyProfile,nextIdx+idxOffset+nPointsPerBodyProfile);
            end
        end
    case 2
        % Quadrilaterals
        for iProfile=1:nBodyProfiles-1
            idxOffset=(iProfile-1)*nPointsPerBodyProfile;
            for iPoint=1:nPointsPerBodyProfile-1
                nextIdx=iPoint+1;
                fprintf(fileID,'%d %d %d %d %d\n',nQuadSides,iPoint+idxOffset,nextIdx+idxOffset,...
                    nextIdx+idxOffset+nPointsPerBodyProfile,iPoint+idxOffset+nPointsPerBodyProfile);
            end
        end    
end
%%%%%%%%%%%%%%%%%%%%%% Pylon polygons %%%%%%%%%%%%%%%%%%%%%%%%
% Triangles (front)
for iPoint=1:nPointsPerPylonProfile-1
    nextIdx=iPoint+1+idxPylonHead;
    fprintf(fileID,'%d %d %d %d\n',nTriangleSides,idxPylonHead,iPoint+idxPylonHead,nextIdx);
end
%  Triangles (rear)
idxOffset=nPointsPerPylonProfile*(nPylonProfiles-1);
for iPoint=1:nPointsPerPylonProfile-1
    nextIdx=iPoint+1+idxPylonHead;
    fprintf(fileID,'%d %d %d %d\n',nTriangleSides,idxPylonTail,iPoint+idxPylonHead+idxOffset,nextIdx+idxOffset);
end
% Pylon polygons (depends on mesh type)
switch(meshType)
    case 1
        % Triangles
        for iProfile=1:nPylonProfiles-1
            idxOffset=(iProfile-1)*nPointsPerPylonProfile;
            for iPoint=1:nPointsPerPylonProfile-1
                nextIdx=iPoint+1+idxPylonHead;
                fprintf(fileID,'%d %d %d %d\n',nTriangleSides,iPoint+idxPylonHead+idxOffset,nextIdx+idxOffset,nextIdx+idxOffset+nPointsPerPylonProfile);
                fprintf(fileID,'%d %d %d %d\n',nTriangleSides,iPoint+idxPylonHead+idxOffset,iPoint+idxPylonHead+idxOffset+nPointsPerPylonProfile,nextIdx+idxOffset+nPointsPerPylonProfile);
            end
        end
    case 2
        % Quadrilaterals
        for iProfile=1:nPylonProfiles-1
            idxOffset=(iProfile-1)*nPointsPerPylonProfile;
            for iPoint=1:nPointsPerPylonProfile-1
                nextIdx=iPoint+1+idxPylonHead;
                fprintf(fileID,'%d %d %d %d %d\n',nQuadSides,iPoint+idxPylonHead+idxOffset,nextIdx+idxOffset,...
                    nextIdx+idxOffset+nPointsPerPylonProfile,iPoint+idxPylonHead+idxOffset+nPointsPerPylonProfile);
            end
        end    
end
% Close file
fclose(fileID);
end

%% Write PLT File
function WritePLT(meshType,bodyHead,bodyTail,x_body,y_body,z_body,pylonHead,pylonTail,x_pylon,y_pylon,z_pylon)
% Specify file name
fileName='RobinFuselagePylon';
% Create and open file
fileID=fopen(strcat(fileName,'.plt'),'w');
fprintf(fileID,'TITLE = "ROBIN Fuselage"\n');
fprintf(fileID,'VARIABLES = "x" "y" "z"\n');

%% WRITE ELEMENTS
%%%%%%%%%%%%%%%%%%%%%%%% Body Elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nProfiles=size(x_body,2);
nPointsPerProfile=size(x_body,1);
% head points (triangles)
fprintf(fileID,'ZONE T="FUSELAGEHEAD" N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n',nPointsPerProfile+1,nPointsPerProfile-1);
fprintf(fileID,'%f %f %f\n',bodyHead(1),bodyHead(2),bodyHead(3));
for iPoint=1:nPointsPerProfile
    fprintf(fileID,'%f %f %f\n',x_body(iPoint,1),y_body(iPoint,1),z_body(iPoint,1));
end
% Indices of nodes of head panels
for iPoint=1:nPointsPerProfile-1
    nextIdx=iPoint+1;
    fprintf(fileID,'%d %d %d\n',1,iPoint+1,nextIdx+1);
end

% % body points (quadrilaterals)
nBodyPoints=nProfiles*nPointsPerProfile;
switch(meshType)
    case 1
        % Triangles
        nBodyPanels=(nProfiles-1)*(nPointsPerProfile-1)*2;
        fprintf(fileID,'ZONE T="FUSELAGEBODY" N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n',nBodyPoints,nBodyPanels);
        for iProfile=1:nProfiles    
            for iPoint=1:nPointsPerProfile
                fprintf(fileID,'%f %f %f\n',x_body(iPoint,iProfile),y_body(iPoint,iProfile),z_body(iPoint,iProfile));
            end
        end
        % Indices of nodes of body panels
        for iProfile=1:nProfiles-1
            idxOffset=(iProfile-1)*nPointsPerProfile;
            for iPoint=1:nPointsPerProfile-1
                nextIdx=iPoint+1;
                fprintf(fileID,'%d %d %d\n',iPoint+idxOffset,nextIdx+idxOffset,nextIdx+idxOffset+nPointsPerProfile);
                fprintf(fileID,'%d %d %d\n',iPoint+idxOffset,nextIdx+idxOffset+nPointsPerProfile,iPoint+idxOffset+nPointsPerProfile);
            end
        end
    case 2
        % Quadrilaterals
        nBodyPanels=(nProfiles-1)*(nPointsPerProfile-1);
        fprintf(fileID,'ZONE T="FUSELAGEBODY" N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n',nBodyPoints,nBodyPanels);
        for iProfile=1:nProfiles    
            for iPoint=1:nPointsPerProfile
                fprintf(fileID,'%f %f %f\n',x_body(iPoint,iProfile),y_body(iPoint,iProfile),z_body(iPoint,iProfile));
            end
        end
        % Indices of nodes of body panels
        for iProfile=1:nProfiles-1
            idxOffset=(iProfile-1)*nPointsPerProfile;
            for iPoint=1:nPointsPerProfile-1
                nextIdx=iPoint+1;
                fprintf(fileID,'%d %d %d %d\n',iPoint+idxOffset,nextIdx+idxOffset,...
                    nextIdx+idxOffset+nPointsPerProfile,iPoint+idxOffset+nPointsPerProfile);
            end
        end
end

% tail points (triangles)
fprintf(fileID,'ZONE T="FUSELAGETAIL" N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n',nPointsPerProfile+1,nPointsPerProfile-1);
for iPoint=1:nPointsPerProfile
    fprintf(fileID,'%f %f %f\n',x_body(iPoint,nProfiles),y_body(iPoint,nProfiles),z_body(iPoint,nProfiles));
end
fprintf(fileID,'%f %f %f\n',bodyTail(1),bodyTail(2),bodyTail(3));
% Indices of nodes of tail panels
for iPoint=1:nPointsPerProfile-1
    nextIdx=iPoint+1;
    fprintf(fileID,'%d %d %d\n',nPointsPerProfile+1,iPoint,nextIdx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Pylon Elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%
nProfiles=size(x_pylon,2);
nPointsPerProfile=size(x_pylon,1);
% head points (triangles)
fprintf(fileID,'ZONE T="FUSELAGEHEAD" N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n',nPointsPerProfile+1,nPointsPerProfile-1);
fprintf(fileID,'%f %f %f\n',pylonHead(1),pylonHead(2),pylonHead(3));
for iPoint=1:nPointsPerProfile
    fprintf(fileID,'%f %f %f\n',x_pylon(iPoint,1),y_pylon(iPoint,1),z_pylon(iPoint,1));
end
% Indices of nodes of head panels
for iPoint=1:nPointsPerProfile-1
    nextIdx=iPoint+1;
    fprintf(fileID,'%d %d %d\n',1,iPoint+1,nextIdx+1);
end

% pylon points (quadrilaterals)
nBodyPoints=nProfiles*nPointsPerProfile;
switch(meshType)
    case 1
        % Triangles
        nBodyPanels=(nProfiles-1)*(nPointsPerProfile-1)*2;
        fprintf(fileID,'ZONE T="FUSELAGEBODY" N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n',nBodyPoints,nBodyPanels);
        for iProfile=1:nProfiles    
            for iPoint=1:nPointsPerProfile
                fprintf(fileID,'%f %f %f\n',x_pylon(iPoint,iProfile),y_pylon(iPoint,iProfile),z_pylon(iPoint,iProfile));
            end
        end
        % Indices of nodes of body panels
        for iProfile=1:nProfiles-1
            idxOffset=(iProfile-1)*nPointsPerProfile;
            for iPoint=1:nPointsPerProfile-1
                nextIdx=iPoint+1;
                fprintf(fileID,'%d %d %d\n',iPoint+idxOffset,nextIdx+idxOffset,nextIdx+idxOffset+nPointsPerProfile);
                fprintf(fileID,'%d %d %d\n',iPoint+idxOffset,nextIdx+idxOffset+nPointsPerProfile,iPoint+idxOffset+nPointsPerProfile);
            end
        end
    case 2
        % Quadrilaterals
        nBodyPanels=(nProfiles-1)*(nPointsPerProfile-1);
        fprintf(fileID,'ZONE T="FUSELAGEBODY" N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n',nBodyPoints,nBodyPanels);
        for iProfile=1:nProfiles    
            for iPoint=1:nPointsPerProfile
                fprintf(fileID,'%f %f %f\n',x_pylon(iPoint,iProfile),y_pylon(iPoint,iProfile),z_pylon(iPoint,iProfile));
            end
        end
        % Indices of nodes of body panels
        for iProfile=1:nProfiles-1
            idxOffset=(iProfile-1)*nPointsPerProfile;
            for iPoint=1:nPointsPerProfile-1
                nextIdx=iPoint+1;
                fprintf(fileID,'%d %d %d %d\n',iPoint+idxOffset,nextIdx+idxOffset,...
                    nextIdx+idxOffset+nPointsPerProfile,iPoint+idxOffset+nPointsPerProfile);
            end
        end
end

% tail points (triangles)
fprintf(fileID,'ZONE T="FUSELAGETAIL" N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n',nPointsPerProfile+1,nPointsPerProfile-1);
for iPoint=1:nPointsPerProfile
    fprintf(fileID,'%f %f %f\n',x_pylon(iPoint,nProfiles),y_pylon(iPoint,nProfiles),z_pylon(iPoint,nProfiles));
end
fprintf(fileID,'%f %f %f\n',pylonTail(1),pylonTail(2),pylonTail(3));
% Indices of nodes of tail panels
for iPoint=1:nPointsPerProfile-1
    nextIdx=iPoint+1;
    fprintf(fileID,'%d %d %d\n',nPointsPerProfile+1,iPoint,nextIdx);
end
% Close file
fclose(fileID);
end