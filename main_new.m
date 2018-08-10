clear; close all; clc;

%% plot settings
fontSize=15;
faceAlpha1=0.3;
faceAlpha2=1;
cMap=gjet(4);
patchColor=cMap(1,:);
markerSize=10;
edgeColor=0.25*ones(1,3);
edgeWidth=1.5;

%% path names
% filePath=mfilename('fullpath');
% savePath=fullfile(fileparts(fileparts(filePath)),'GIBBON-master','data','temp');
% modelName=fullfile(savePath,'sphere_new');
% modelName=fullfile(savePath,'sphere3');

% build a sphere surface
% r1=1; %Outer sphere radius
% numRefine=2; %Number of refinement steps from icosahedron
% faceBoundMarker=3; %Face marker for outer sphere
% 
% [Fq,Vq,~]=geoSphere(numRefine,r1);
% 
% faceBoundaryMarker_q=faceBoundMarker*ones(size(Fq,1),1); %Create boundary markers for faces
% 
% hf=cFigure;
% title('Surface models','FontSize',fontSize);
% xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
% hold on;
% 
% patch('Faces',Fq,'Vertices',Vq,'FaceColor','flat','CData',faceBoundaryMarker_q,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
% [hp]=patchNormPlot(Fq,Vq,0.25);
% 
% colormap(autumn(2));
% colorbar;
% camlight headlight;
% set(gca,'FontSize',fontSize);
% view(3); axis tight;  axis equal;  grid on;
% 
% %CREATING A SOLID TETRAHEDRAL MESH USING TETGEN
% % %% CREATING THE INPUT STRUCTURE
% [regionA]=tetVolMeanEst(Fq,Vq); %Volume for regular tets
% 
% stringOpt='-pq1.2AaYQ';
% 
% inputStruct.stringOpt=stringOpt;
% inputStruct.Faces=Fq;
% inputStruct.Nodes=Vq;
% inputStruct.holePoints=[];
% inputStruct.faceBoundaryMarker=faceBoundaryMarker_q; %Face boundary markers
% inputStruct.regionPoints=[0 0 0]; %region points
% inputStruct.regionA=regionA;
% inputStruct.minRegionMarker=2; %Minimum region marker
% inputStruct.modelName=modelName;
% 
% %% Mesh model using tetrahedral elements using tetGen
% [meshOutput]=runTetGen(inputStruct); %Run tetGen
% 
% Fb=meshOutput.facesBoundary;
% Cb=meshOutput.boundaryMarker;
% V=meshOutput.nodes;
% CE=meshOutput.elementMaterialID;
% El=meshOutput.elements;

% meshView(meshOutput,[]);

% defaultFolder = fileparts(fileparts(mfilename('fullpath')));
% pathName=fullfile(defaultFolder,'GibbonCode','data','STL');
% fileName=fullfile(pathName,'brain_decimated.stl');
[A] = textread('/home/x17wang/Bureau/MATLAB/GibbonCode/data/MESH/sphere7.mesh');

V = A(2:A(1,1)+1,1:3);
El = A(A(1,1)+3:A(A(1,1)+2,1)+A(1,1)+2,2:5);
Fb = A(A(A(1,1)+2,1)+A(1,1)+4:end,2:4);

% [stlStruct] = import_STL(fileName);
% 
% F=stlStruct.solidFaces{1};
% V=stlStruct.solidVertices{1};

% find center of mass and dimension of the mesh
maxx = -1e9;
minx = 1e9;
maxy = -1e9;
miny = 1e9;
maxz = -1e9;
minz = 1e9;
cog = zeros(1,3);
for i = 1:size(V,1)
    maxx = max(maxx, V(i,1)); minx = min(minx,V(i,1));
    maxy = max(maxy, V(i,2)); miny = min(miny,V(i,2));
    maxz = max(maxz, V(i,3)); minz = min(minz,V(i,3));
    cog = cog + V(i,:);
end;
cog = cog/size(V,1);
maxd = max(max(max(abs(maxx-cog(1)), abs(minx-cog(1))), max(abs(maxy-cog(2)), abs(miny-cog(2)))), max(abs(maxz-cog(3)), abs(minz-cog(3))));

% change mesh information by values normalized 
for i = 1:size(V,1)
    V(i,1) = (V(i,1) - cog(1))/maxd;
    V(i,2) = (V(i,2) - cog(2))/maxd;
    V(i,3) = -(V(i,3) - cog(3))/maxd;
end;

% selecting half of the model to see interior
% Y=V(:,2); YE=mean(Y(El),2);
% logicCutView=YE>mean(Y);
% [Fs,Cs]=element2patch(El(logicCutView,:),CE(logicCutView),'tet4');

% cut view of tetrahedral mesh model
% cFigure;
% hold on;
% title('Cut view of tetrahedral mesh model','FontSize',fontSize);
% gpatch(Fb,V,0.5*ones(1,3),'none',faceAlpha1);
% gpatch(Fs,V,Cs,'k',faceAlpha2);
% plotV(V(unique(Fs(:)),:),'k.','MarkerSize',markerSize);
% camlight headlight;
% axisGeom(gca,fontSize);
% axis off;Oui
% colormap(autumn);
% drawnow;

% the outer surface normals
% hf=cFigure;
% title('The outer surface normals','FontSize',fontSize);
% xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
% hold on;
% patch('Faces',Fb,'Vertices',V,'FaceColor','flat','CData',Cb,'FaceAlpha',0.5);
% [hp]=patchNormPlot(Fb,V,0.5); 
% colormap jet; colorbar;
% camlight headlight;
% set(gca,'FontSize',fontSize);
% view(3); axis tight;  axis equal;  grid on;

% % Merging nodes (nodes are not merged in stl)
% [~,ind1,ind2]=unique(pround(V,5),'rows');
% V=V(ind1,:);
% F=ind2(F);
% 
cFigure; hold on;
title('Surface model','FontSize',fontSize);
gpatch(Fb,V,patchColor,'k',faceAlpha1);
camlight headlight;
axisGeom(gca,fontSize);
drawnow;
% 
% faceBoundaryMarker=ones(size(F,1),1);
% 
% [V_regions]=getInnerPoint(F,V);
% 
% V_holes=[];
% 
% [regionA]=tetVolMeanEst(F,V); %Volume for regular tets
% 
% stringOpt='-pq1.2AaYQ';
% 
% inputStruct.stringOpt=stringOpt;
% inputStruct.Faces=F;
% inputStruct.Nodes=V;
% inputStruct.holePoints=V_holes;
% inputStruct.faceBoundaryMarker=faceBoundaryMarker; %Face boundary markers
% inputStruct.regionPoints=V_regions; %region points
% inputStruct.regionA=regionA;
% inputStruct.minRegionMarker=2; %Minimum region marker
% inputStruct.modelName=modelName;
% 
% [meshOutput]=runTetGen(inputStruct); %Run tetGen
% 
% Fb=meshOutput.facesBoundary;
% Cb=meshOutput.boundaryMarker;
% V=meshOutput.nodes;
% CE=meshOutput.elementMaterialID;
% El=meshOutput.elements;
% 
% meshView(meshOutput,[]);

%% CONTROL PARAMETERS
% FEA control settings
% numTimeSteps=20; %Number of time steps desired
% max_refs=25; %Max reforms
% max_ups=0; %Set to zero to use full-Newton iterations per time step
% opt_iter=10; %Optimum number of iterations
% max_retries=25; %Maximum number of retries
% dtmin=(1/numTimeSteps)/100; %Minimum time step size
% dtmax=1/numTimeSteps; %Maximum time step size
% t_load=0.5; %Time from start to max load
% t_unload=0.5;  %Time from max load to end
% t_wait=1; %Additional wait time
% t_total=t_load+t_unload+t_wait; %Total simulation time
% t_step_ini=0.05; %Initial desired step size
% t_step_max=t_step_ini; %Maximum step size
% numTimeSteps=round(t_total/t_step_ini);
% t_step=t_total/numTimeSteps;
% 
% uncoupledLaw=1; %1=uncoupled, 2=coupled
% 
% numTimeSteps=20; %Number of time steps desired
% max_refs=25; %Max reforms
% max_ups=0; %Set to zero to use full-Newton iterations
% opt_iter=10; %Optimum number of iterations
% max_retries=5; %Maximum number of retires
% dtmin=(1/numTimeSteps)/100; %Minimum time step size
% dtmax=1/numTimeSteps; %Maximum time step size
% 
% r1=0.52; %Outer sphere radius
% numRefine=3; %Number of refinement steps from icosahedron
% faceBoundMarker=3; %Face marker for outer sphere
% 
% [Fq,Vq,~]=geoSphere(numRefine,r1);
% 
% faceBoundaryMarker_q=faceBoundMarker*ones(size(Fq,1),1); %Create boundary markers for faces
% 
% hf=cFigure;
% title('Surface models','FontSize',fontSize);
% xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
% hold on;
% 
% patch('Faces',Fq,'Vertices',Vq,'FaceColor','flat','CData',faceBoundaryMarker_q,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
% [hp]=patchNormPlot(Fq,Vq,0.25);
% 
% colormap(autumn(2));
% colorbar;
% camlight headlight;
% set(gca,'FontSize',fontSize);
% view(3); axis tight;  axis equal;  grid on;
% 
% %CREATING A SOLID TETRAHEDRAL MESH USING TETGEN
% % %% CREATING THE INPUT STRUCTURE
% [regionA]=tetVolMeanEst(Fq,Vq); %Volume for regular tets
% 
% stringOpt='-pq1.2AaYQ';
% 
% inputStruct.stringOpt=stringOpt;
% inputStruct.Faces=Fq;
% inputStruct.Nodes=Vq;
% inputStruct.holePoints=[];
% inputStruct.faceBoundaryMarker=faceBoundaryMarker_q; %Face boundary markers
% inputStruct.regionPoints=[0 0 0]; %region points
% inputStruct.regionA=regionA;
% inputStruct.minRegionMarker=2; %Minimum region marker
% inputStruct.modelName=modelName;
% 
% %% Mesh model using tetrahedral elements using tetGen
% [meshOutput]=runTetGen(inputStruct); %Run tetGen
% 
% Fb=meshOutput.facesBoundary;
% Cb=meshOutput.boundaryMarker;
% V=meshOutput.nodes;
% CE=meshOutput.elementMaterialID;
% El=meshOutput.elements;
% 
% % Selecting half of the model to see interior
% Y=V(:,2); YE=mean(Y(El),2);
% logicCutView=YE>mean(Y);
% [Fs,Cs]=element2patch(El(logicCutView,:),CE(logicCutView),'tet4');
% 
% % Cut view of tetrahedral mesh model
% cFigure;
% hold on;
% title('Cut view of tetrahedral mesh model','FontSize',fontSize);
% gpatch(Fb,V,0.5*ones(1,3),'none',faceAlpha1);
% gpatch(Fs,V,Cs,'k',faceAlpha2);
% plotV(V(unique(Fs(:)),:),'k.','MarkerSize',markerSize);
% camlight headlight;
% axisGeom(gca,fontSize);
% axis off;
% colormap(autumn);
% drawnow;
% 
% % The outer surface normals
% hf=cFigure;
% title('The outer surface normals','FontSize',fontSize);
% xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
% hold on;
% patch('Faces',Fb,'Vertices',V,'FaceColor','flat','CData',Cb,'FaceAlpha',0.5);
% [hp]=patchNormPlot(Fb,V,0.5); 
% colormap jet; colorbar;
% camlight headlight;
% set(gca,'FontSize',fontSize);
% view(3); axis tight;  axis equal;  grid on;

% %Define BC's
% forceMagnitude=[1 1 1];
% bcPrescribeList = unique(Fb(:));
% bcPrescribeMagnitude = forceMagnitude(ones(1,numel(bcPrescribeList)),:);

% Build a quadrilateral surface
% boxDim=[2 5 2];  %Dimensions
% boxEl=[2 3 2]; %Number of elements
% [Fq,Vq,faceBoundaryMarker_q]=quadBox(boxDim,boxEl);
% 
% [V_regions]=getInnerPoint(Fq,Vq); % Define region points
% V_holes=[]; % Define hole points
% 
% % PlottVing surface models
% hf=cFigure;
% title('Surface models','FontSize',fontSize);
% xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
% hold on;
% patch('Faces',Fq,'Vertices',Vq,'FaceColor','flat','CData',faceBoundaryMarker_q,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
% colormap(autumn(2));
% colorbar;
% camlight headlight;
% set(gca,'FontSize',fontSize);
% view(3); axis tight;  axis equal;  grid on;
% 
% %% CREATING THE INPUT STRUCTURE
% [regionA]=tetVolMeanEst(Fq,Vq); %Volume for regular tets
% 
% stringOpt='-pq1.2AaYQ';
% 
% inputStruct.stringOpt=stringOpt;
% inputStruct.Faces=Fq;
% inputStruct.Nodes=Vq;
% inputStruct.holePoints=[];
% inputStruct.faceBoundaryMarker=faceBoundaryMarker_q; %Face boundary markers
% % inputStruct.regionPoints=[0 0 0]; %region points
% inputStruct.regionPoints=V_regions; %region points
% inputStruct.regionA=regionA;
% inputStruct.minRegionMarker=2; %Minimum region marker
% inputStruct.modelName=modelName;
% 
% %% Mesh model using tetrahedral elements using tetGen
% [meshOutput]=runTetGen(inputStruct); %Run tetGen
% 
% Fb=meshOutput.facesBoundary;
% Cb=meshOutput.boundaryMarker;
% V=meshOutput.nodes;
% CE=meshOutput.elementMaterialID;
% El=meshOutput.elements;
% 
% % Selecting half of the model to see interior
% Y=V(:,2); YE=mean(Y(El),2);
% logicCutView=YE>mean(Y);
% [Fs,Cs]=element2patch(El(logicCutView,:),CE(logicCutView),'tet4');
% 
% % Cut view of tetrahedral mesh model
% cFigure;
% hold on;
% title('Cut view of tetrahedral mesh model','FontSize',fontSize);
% gpatch(Fb,V,0.5*ones(1,3),'none',faceAlpha1);
% gpatch(Fs,V,Cs,'k',faceAlpha2);
% plotV(V(unique(Fs(:)),:),'k.','MarkerSize',markerSize);
% camlight headlight;
% axisGeom(gca,fontSize);
% axis off;
% colormap(autumn);
% drawnow;
%  
% % The outer surface normals
% hf=cFigure;
% title('The outer surface normals','FontSize',fontSize);
% xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
% hold on;
% patch('Faces',Fb,'Vertices',V,'FaceColor','flat','CData',Cb,'FaceAlpha',0.5);
% [hp]=patchNormPlot(Fb,V,0.5); 
% colormap jet; colorbar;
% camlight headlight;
% set(gca,'FontSize',fontSize);
% view(3); axis tight;  axis equal;  grid on;

% Mark non-growing areas
V0=V;
gr=zeros(1,size(V,1));
parfor i=1:size(V,1)
    qp = V0(i,:);
    X = [(qp(1)+0.1)*0.714, qp(2), qp(3)-0.05];
    rqp = length(X);
    if rqp < 0.6
        gr(i) = max(1.0-10.0*(0.6-rqp),0.0);
    else
        gr(i) = 1.0;
    end
end    

%% Determine surface nodes and index maps
ns(size(V,1))=0;
parfor j=1:size(V,1)
    for i=1:size(Fb,1)
        for k=1:size(Fb,2)
            if j == Fb(i,k) 
                ns(j)=1;
            end
        end
    end
end

snb(size(V,1))=0;p=1;sn=0;
for i=1:size(V,1)
   if ns(i)==1
        sn(p)=i;   % surface to full mesh
        snb(i)=p;  % full mesh to surface
        p=p+1;
   end
end

% Normal vector of the surface points
no=cell(1,size(sn,2)); % all the points of surface
n1=cell(1,size(sn,2));
for i=1:size(Fb,1)
    for j=1:size(Fb,2)
        no{snb(Fb(i,j))}= [0 0 0];
    end
end

% No=cell(1,size(sn,2));
for i=1:size(Fb,1)
%    for j=1:size(Fb,2)
   No = cross(V0(Fb(i,2),:)-V0(Fb(i,1),:), V0(Fb(i,3),:)-V0(Fb(i,1),:));
   no{snb(Fb(i,1))}= no{snb(Fb(i,1))}+No;       % Normal for each surface point
   no{snb(Fb(i,2))}= no{snb(Fb(i,2))}+No; 
   no{snb(Fb(i,3))}= no{snb(Fb(i,3))}+No; 
   %        n1{snb(Fb(i,j))}= no{snb(Fb(i,j))}./norm(no{snb(Fb(i,j))});
%    end                                                           %Attention: index is surface index = snb(full mesh index)
end

parfor i=1:size(sn,2)
    n1{i}= no{i}/norm(no{i});
end

nsn=length(sn); % number of surface nodes
% Find nearst point
csn=zeros(1,size(V0,1));
for i=1:size(V0,1)
    d2min=1e9;
    for j=1:nsn     %all the surface points,j:surface index B{i}(isnan(B{i}) | isinf(B{i}))=1;
        d2=dot((V0(sn(j),:)-V0(i,:)),(V0(sn(j),:)-V0(i,:)));        %calcul distance 
        if d2<d2min
            d2min=d2;
            q=j;
        end
    end
    csn(i)=q;
    d2s(i)=sqrt(d2min);
%     else
%         csn(i)=snb(i);      % the nearst point is itself
%         d2s(i)=0;           % distance = 0
%     end    
end


% Find the normals for each tetrahedron
NL=cell(size(El,1),4);
NL_TOTAL=cell(size(El,1),1);
for i = 1:size(El,1)
    NE(i,1) = csn(El(i,1));
    NE(i,2) = csn(El(i,2));
    NE(i,3) = csn(El(i,3));
    NE(i,4) = csn(El(i,4));
    NL{i,1} = n1{NE(i,1)};
    NL{i,2} = n1{NE(i,2)};
    NL{i,3} = n1{NE(i,3)};
    NL{i,4} = n1{NE(i,4)};
    NL_TOTAL{i,1} = NL{i,1}+NL{i,2}+NL{i,3}+NL{i,4};
    NL_TOTAL{i,1} = NL_TOTAL{i,1}/norm(NL_TOTAL{i,1});
end

% Define BC's
bcPrescribeList=[1:size(V,1)]';

%% loop for iterations of matlab and steps of FEBio
% p=1;  % total time of deformation
n=12; %number of iterations for calculating the forces
t=-0.25; %current time
dt=0.050*sqrt(0.0025*0.01*0.01/5.0); %time step size for matlab
dtt=60*0.050*sqrt(0.0025*0.01*0.01/5.0); %time step size for FEBio
% dt=1/n;
% n=round(1/dt); % number of time step 
% dt=1/n; %time step size
d=1; %for creating the files.feb,.log,.xplt
Vt=zeros(3,size(V,1));
Vtold=zeros(size(unique(Fb),1),3);
NNLt=cell(1,nsn);
ub=0; vb=0; wb=0;

G=cell(1,size(El,1));
for i = 1:size(El,1)
    G{i}= [1 0 0;0 1 0;0 0 1];
end;
% E=cell(1,n);
% v=cell(1,n);
GG=cell(1,n);
Ftt=cell(1,n);
Fc=cell(1,n);
for i = 1:n
    [W,E,v,Ft,Ftt{i},G]=definition_material(V,V0,El,Fb,dt,Vtold,G,t,Vt,sn,gr,NNLt,ub,vb,wb,d2s,NL_TOTAL,nsn);
    Fc{i}=Ft;
    GG{i}=G;
    % Calcul of the positions of centre-of-gravity of tetrahedrons for plotting
% the normals of tetrahedrons
%     for i = 1:size(El,1)
%         C(i,1) = (V(El(i,1),1)+V(El(i,2),1)+V(El(i,3),1)+V(El(i,4),1))/4.0;
%         C(i,2) = (V(El(i,1),2)+V(El(i,2),2)+V(El(i,3),2)+V(El(i,4),2))/4.0;
%         C(i,3) = (V(El(i,1),3)+V(El(i,2),3)+V(El(i,3),3)+V(El(i,4),3))/4.0;
%     end;
%     hf=cFigure;
%     title('Growth Tensors','FontSize',fontSize);
%     xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
%     % tetramesh(El,V);
%     hold on;
%     for i = 1:size(El,1) 
%      NNLt   quiver3(C(i,1),C(i,2),C(i,3),G{i}(1,1),G{i}(1,2),G{i}(1,3),2,'linewidth',4);
%         quiver3(C(i,1),C(i,2),C(i,3),G{i}(2,1),G{i}(2,2),G{i}(2,3),2,'linewidth',4);
%         quiver3(C(i,1),C(i,2),C(i,3),G{i}(3,1),G{i}(3,2),G{i}(3,3),2,'linewidth',4);
%     end;
%     % hp=patch('Faces',El,'Vertices',V,'FaceColor','flat','CData',faceBoundaryMarker_q,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
%     % set(hp,'EdgeColor','none','FaceColor','k');
%     patch('Faces',Fb,'Vertices',V,'FaceColor','flat','CData',Cb,'FaceAlpha',0.5);
% %     [hp]=patchNormPlot(Fb,V,0.5);
%     colormap(autumn(2));
%     colorbar;
%      headlight;
%     set(gca,'FontSize',fontSize);
%     view(3); axis tight;  axis equal;  grid on;% a = [];
    %Define BC's
    parfor j = 1:size(V,1)
%         forceMagnitude=[1 1 1];
        forceMagnitude=[Ft(1,j) Ft(2,j) Ft(3,j)];
%         bcPrescribeList = unique(Fb(:));
%     bcPrescribeMagnitude = forceMagnitude(ones(1,numel(bcPrescribeList)),:);
        bcPrescribeMagnitude(j,:)=forceMagnitude(:);
    end
    [V]=runFEBio_new(V,El,E,v,W,d,bcPrescribeList,bcPrescribeMagnitude,dtt);
    d=d+1;
end