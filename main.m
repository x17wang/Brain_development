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

%% Path names
filePath=mfilename('fullpath');
savePath=fullfile(fileparts(fileparts(filePath)),'GibbonCode','data','temp');
% modelName=fullfile(savePath,'sphere_pos');
modelName=fullfile(savePath,'sphere2');
% 
% defaultFolder = fileparts(fileparts(mfilename('fullpath')));
% pathName=fullfile(defaultFolder,'GibbonCode','data','STL');
% fileName=fullfile(pathName,'brain_decimated.stl');
[A] = textread('/home/xiaoyu/Desktop/MATLAB/GibbonCode/data/MESH/sphere2.mesh');
% V = A(2:84,1:3);
% El = A(86:266,2:5);
% Fb = A(268:end,2:4);
V = A(2:35590,1:3);
El = A(35592:229127,2:5);
Fb = A(229129:end,2:4);
% [stlStruct] = import_STL(fileName);
% 
% F=stlStruct.solidFaces{1};
% V=stlStruct.solidVertices{1};
% 
% 
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
% % Plotting surface models
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

%Define BC's
% forceMagnitude=[1 1 1];

%% loop for iterations of matlab and steps of FEBio
p=1;  % total time of deformation
n=20; % number of time step 0
% dt = 0.050*sqrt(0.0025*0.01*0.01/5.0); %time step size
dt=p/n; %time step size
t=1; %for creating the files.feb,.log,.xplt
V0=V;
Vtold=zeros(size(unique(Fb),1),3);
P=cell(1,n);
G = cell(1,size(El,1));
for i = 1:size(El,1)
    G{i}= [1 0 0;0 1 0;0 0 1];
end;
% E=cell(1,n);
% v=cell(1,n);
for i = 1:n
    [W,E,v,Ft]=definition_material(V,V0,El,Fb,dt,i,Vtold,G);
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
%         quiver3(C(i,1),C(i,2),C(i,3),G{i}(1,1),G{i}(1,2),G{i}(1,3),2,'linewidth',4);
%         quiver3(C(i,1),C(i,2),C(i,3),G{i}(2,1),G{i}(2,2),G{i}(2,3),2,'linewidth',4);
%         quiver3(C(i,1),C(i,2),C(i,3),G{i}(3,1),G{i}(3,2),G{i}(3,3),2,'linewidth',4);
%     end;
%     % hp=patch('Faces',El,'Vertices',V,'FaceColor','flat','CData',faceBoundaryMarker_q,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
%     % set(hp,'EdgeColor','none','FaceColor','k');
%     patch('Faces',Fb,'Vertices',V,'FaceColor','flat','CData',Cb,'FaceAlpha',0.5);
% %     [hp]=patchNormPlot(Fb,V,0.5);
%     colormap(autumn(2));
%     colorbar;
%     camlight headlight;
%     set(gca,'FontSize',fontSize);
%     view(3); axis tight;  axis equal;  grid on;% a = [];
    %Define BC's
    for j = 1:size(V,1)
%         forceMagnitude=[1 1 1];
        forceMagnitude=[Ft(1,j) Ft(2,j) Ft(3,j)];
%         bcPrescribeList = unique(Fb(:));
%     bcPrescribeMagnitude = forceMagnitude(ones(1,numel(bcPrescribeList)),:);
        bcPrescribeList=[1:size(V,1)]';
        bcPrescribeMagnitude(j,:)=forceMagnitude(:);
    end
    [V]=runFEBio(V,El,E,v,W,t,modelName,bcPrescribeList,bcPrescribeMagnitude);
    t=t+1;
end