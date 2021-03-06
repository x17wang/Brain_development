function [W,E,v,Ft,Ftt,G]=definition_material(V,V0,El,Fb,dt,Vtold,G,t,Vt,sn,gr,NNLt,ub,vb,wb,d2s,NL_TOTAL,nsn)

%% Determine at
% t = -0.35+dt*n; % time to calculate at
gamma = 0.5;
rho = 0.0025;
at = 1.7+7.4*t;% growth tensor changes with at changes as time goes on
H = 0.031+0.01*t; % Un deformed thickness of the gray matter
% if t>=0.0
%     at = 1.85-1.85*t;
%     H = 0.042+0.042*t;
% end;
% H = 0.041 + 0.01*t;


% Calcul of the positions of centre-of-gravity of tetrahedrons for plotting
% the normals of tetrahedrons
% for i = 1:size(El,1)
%     C(i,1) = (V(El(i,1),1)+V(El(i,2),1)+V(El(i,3),1)+V(El(i,4),1))/4.0;
%     C(i,2) = (V(El(i,1),2)+V(El(i,2),2)+V(El(i,3),2)+V(El(i,4),2))/4.0;
%     C(i,3) = (V(El(i,1),3)+V(El(i,2),3)+V(El(i,3),3)+V(El(i,4),3))/4.0;
% end;

%%
%         hf=cFigure;
%         title('The tetrahedral normals','FontSize',fontSize);
%         xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
%         % tetramesh(El,V);
%         hold on;
%         Ns = cell2mat(NL_TOTAL);
%         quiver3(C(:,1),C(:,2),C(:,3),Ns(:,1),Ns(:,2),Ns(:,3),1,'linewidth',5);
%         hp=patch('Faces',El,'Vertices',V,'FaceColor','flat','CData',faceBoundaryMarker_q,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
%         set(hp,'EdgeColor','none','FaceColor','k');
%         patch('Faces',Fb,'Vertices',V,'FaceColor','flat','CData',Cb,'FaceAlpha',0.5);
%         [hp]=patchNormPlot(Fb,V,0.5);
%         colormap(autumn(2));
%         colorbar;
%         camlight headlight;
%         set(gca,'FontSize',fontSize);
%         view(3); axis tight;  axis equal;  grid on;% a = [];
% [Fni,Vni,~]=quiver3Dpatch(El(:,1),(:,2),V(:,3),Ns(:,1),Ns(:,2),Ns(:,3),[],a); %Derive quiver patch data

%% Calculate the contact force
maxDist = 0.0;
% a = 5; %Mesh spacing - set manually based on the average spacing in the mesh
% bw = 200; %Width of a bounding box, centered at origin, that encloses the whole geometry even after growth
% K = 5.0; % bulk modulus
% mw = 2*a; 
a = 0.01; %Mesh spacing - set manually based on the average spacing in the mesh
bw = 3.2; %Width of a bounding box, centered at origin, that encloses the whole geometry even after growth
K = 5.0; % bulk modulus
mw = 8.0*a;
% mw = a; %Width of a cell in the linked cell algorithm for proximity detection
% hs = 0.8*a; 
% hc = 0.5*a;
hs = 0.6*a; %Thickness of proximity skin
hc = 0.2*a; %Thickness of repulsive skin
kc = 10.0*K; %Contact stiffness
nf = size(Fb,1);
parfor i = 1:nsn
    maxDist = max(maxDist, length(V(sn(i),:) - Vtold(i,:)));
end;
if maxDist > 0.5*(hs-hc)
    [NNLt] = createNNLtriangle(NNLt, V', Fb', sn, nsn, nf, hs, bw, mw);
    for i = 1:nsn
        Vtold(i,:) = V(sn(i),:);
    end;
end;

Ft=zeros(3,size(V,1));

for i = 1:nsn
%     if size(NNLt{i},2) >= 1
    for tp = 1:size(NNLt{i},2)
        pt = sn(i);
        tri = NNLt{i}(tp)+1;
        [cc, ub, vb, wb] = closestPointTriangle(V(pt,:), V(Fb(tri,1),:), V(Fb(tri,2),:), V(Fb(tri,3),:), ub, vb, wb);
        cc = cc - V(pt,:);
        rc = length(cc);
        if (rc < hc) && (gr(pt) + gr(Fb(tri,1)) > 0.0)
            cc = cc/norm(cc);
            Ntri = cross(V(Fb(tri,2),:) - V(Fb(tri,1),:), V(Fb(tri,3),:) - V(Fb(tri,1),:));
            Ntri = Ntri/norm(Ntri);
            fn = cc*(rc-hc)/hc*kc*a*a;
            if dot(fn, Ntri) < 0.0
                fn = fn-Ntri*dot(fn,Ntri)*2.0;
            end;
            Ft(:,Fb(tri,1)) = Ft(:,Fb(tri,1))-(fn*ub)';
            Ft(:,Fb(tri,2)) = Ft(:,Fb(tri,2))-(fn*vb)';
            Ft(:,Fb(tri,3)) = Ft(:,Fb(tri,3))-(fn*wb)';
            Ft(:,pt) = Ft(:,pt)+fn';
        end;
    end;
%     end;
end;

Ftt = Ft;

%% Volume of initial tetraedron and deformed tetraedron
Ar = cell(1,size(El,1));
A = cell(1,size(El,1));
V_d = V; 
parfor i = 1:size(V,1)
    Vn0(i) = 0.0;
    Vn(i) = 0.0;
end;

for i = 1:size(El,1)
    % Undeformed
    xr1 = V0(El(i,2),:)-V0(El(i,1),:);
    xr2 = V0(El(i,3),:)-V0(El(i,1),:);
    xr3 = V0(El(i,4),:)-V0(El(i,1),:);
    Ar{i} = [xr1;xr2;xr3]';
    Ar{i} = G{i}*Ar{i};
    
    vol0 = det(Ar{i})/6.0;
    Vn0(El(i,1))= Vn0(El(i,1)) + vol0/4.0;
    Vn0(El(i,2))= Vn0(El(i,2)) + vol0/4.0;
    Vn0(El(i,3))= Vn0(El(i,3)) + vol0/4.0;
    Vn0(El(i,4))= Vn0(El(i,4)) + vol0/4.0;
    
    % Deformed
    x1 = V_d(El(i,2),:)-V_d(El(i,1),:);
    x2 = V_d(El(i,3),:)-V_d(El(i,1),:);
    x3 = V_d(El(i,4),:)-V_d(El(i,1),:);
    A{i} = [x1;x2;x3]'; 
    
    vol = det(A{i})/6.0;
    Vn(El(i,1))= Vn(El(i,1)) + vol/4.0;
    Vn(El(i,2))= Vn(El(i,2)) + vol/4.0;
    Vn(El(i,3))= Vn(El(i,3)) + vol/4.0;
    Vn(El(i,4))= Vn(El(i,4)) + vol/4.0;
end;

%% Deformation 
mug = 1.0; % Shear modulus of gray matter
muw = 1.167; % Shear modulus of white matter
N1 = cell(1,size(El,1));
N2 = cell(1,size(El,1));
N3 = cell(1,size(El,1));
N4 = cell(1,size(El,1));
F = cell(1,size(El,1));
J = cell(1,size(El,1));
B = cell(1,size(El,1));
k = 0.0;
T = cell(1,size(El,1));
P = cell(1,size(El,1));
Ns = cell2mat(NL_TOTAL);
for i = 1:size(El,1)
    % Determine tangential growth profile gm and shear modulus
    gm(i) = 1.0/(1.0 + exp(10.0*(0.25*(d2s(El(i,1))+d2s(El(i,2))+d2s(El(i,3))+d2s(El(i,4)))/H -1.0)))*0.25*(gr(El(i,1))+gr(El(i,2))+gr(El(i,3))+gr(El(i,4)));
    wm(i) = 1.0 - gm(i);
    mu(i) = muw*wm(i) + mug*gm(i);
    
    % Basis vector of reference state
    xr1 = V0(El(i,2),:)-V0(El(i,1),:);
    xr2 = V0(El(i,3),:)-V0(El(i,1),:);
    xr3 = V0(El(i,4),:)-V0(El(i,1),:);
    Ar{i} = [xr1;xr2;xr3]'; % Reference state
    Ar{i} = G{i}*Ar{i}; %Apply growth to reference state
    
    % Undeformed normals
    xr1 = [Ar{i}(1,1) Ar{i}(2,1) Ar{i}(3,1)];
    xr2 = [Ar{i}(1,2) Ar{i}(2,2) Ar{i}(3,2)];
    xr3 = [Ar{i}(1,3) Ar{i}(2,3) Ar{i}(3,3)];
    N1{i} = cross(xr3, xr1);
    N2{i} = cross(xr2, xr3);
    N3{i} = cross(xr1, xr2);
    N4{i} = cross(xr2-xr3, xr1-xr3);
    
    % Deformed basis vectors
    x1 = V_d(El(i,2),:)-V_d(El(i,1),:);
    x2 = V_d(El(i,3),:)-V_d(El(i,1),:);
    x3 = V_d(El(i,4),:)-V_d(El(i,1),:);
    A{i} = [x1;x2;x3]';
    
    F{i} = A{i}*(Ar{i})^(-1); % Deformation gradient
    B{i} = F{i}*(F{i})'; %Left Cauchy-Green strain tensor
    J{i} = det(F{i}); % Relative volume change
    J1(i) = Vn(El(i,1))/Vn0(El(i,1));
    J2(i) = Vn(El(i,2))/Vn0(El(i,2));
    J3(i) = Vn(El(i,3))/Vn0(El(i,3));
    J4(i) = Vn(El(i,4))/Vn0(El(i,4));
    Ja(i) = (J1(i)+J2(i)+J3(i)+J4(i))/4.0; % Averaged nodal volume change
    
    % Defining Young's modulus and Poisson's ratio
    E(i) = (9*K*mu(i))/(3*K+mu(i));
    v(i) = ((3*K)-(2*mu(i)))/((6*K)+(2*mu(i)));
    lam(i) = K-(2*mu(i))/3;
    
    % Calculate the volumetric strain energy density and elastic force
    [ll1,ll2,ll3] = EV(B{i});
    if (ll3>=eps*eps) && (J{i}>0.0) % No need for SVD
%     if J{i}>0.0
        T{i} = (B{i}-eye(3)*trace(B{i})/3.0)*mu(i)/(J{i}*(J{i}^(2.0/3.0))) + eye(3)*K*(Ja(i)-1.0);
        P{i} = T{i}*(inv((F{i})'))*J{i}; %fo,r calculating the forces
%         W(i) = 0.5*mu(i)*(trace(B{i})/J{i}^(2.0/3.0)-3.0)+0.5*K*((J1(i)-1.0)*(J1(i)-1.0)+(J2(i)-1.0)*(J2(i)-1.0)+(J3(i)-1.0)*(J3(i)-1.0)+(J4(i)-1.0)*(J4(i)-1.0))*0.25;
        W_T(i) = 0.5*E(i)/(2*(1+v(i)))*(trace(B{i})/J{i}^(2.0/3.0)-3.0)+0.5*E(i)/(3*(1-2*v(i)))*((J1(i)-1.0)*(J1(i)-1.0)+(J2(i)-1.0)*(J2(i)-1.0)+(J3(i)-1.0)*(J3(i)-1.0)+(J4(i)-1.0)*(J4(i)-1.0))*0.25;
%         W(i) = mu(i)*(((trace(B{i}) - 3)/2) - log(J{i})) + lam(i)*((log(J{i}))^2)/2;
        W(i) = rho;
    else % need SVD
        D = (F{i})'*F{i};
%         D(isnan(D) | isinf(D))=1;
        [U,S] = eig(D);
%         d = zeros(1,3);
%         [U,S] = Eigensystem(D,d);
        l1 = sqrt(S(1,1));
        l2 = sqrt(S(2,2));
        l3 = sqrt(S(3,3));
        
        if det(U)<0.0
            U(1,1) = -U(1,1);
            U(2,1) = -U(2,1);
            U(3,1) = -U(3,1);
        end;
        
        Fdi=zeros(size(D));
        if l1>=10^(-25)
            Fdi(1,1) = 1.0/l1;
            Fdi(2,2) = 1.0/l2;
            Fdi(3,3) = 1.0/l3;
        end;
        
        U1 = F{i}*(U*Fdi);
            
        if l1<10^(-25)
           U1(1,1) = U1(2,2)*U1(3,3)-U1(3,2)*U1(2,3);
           U1(2,1) = U1(3,2)*U1(1,3)-U1(1,2)*U1(3,3);
           U1(3,1) = U1(1,2)*U1(2,3)-U1(2,2)*U1(1,3);
        end;  
        if det(F{i})<0.0
            l1 = -l1;
            U1(1,1) = -U1(1,1);
            U1(2,1) = -U1(2,1);
            U1(3,1) = -U1(3,1);
        end;
        
        Pd=zeros(size(D));
        pow23 = (eps*l2*l3)^(2.0/3.0);
        Pd(1,1) = mu(i)/3.0*(2.0*eps - l2*l2/eps - l3*l3/eps)/pow23 + k*(l1-eps) + K*(Ja(i)-1.0)*l2*l3;
        Pd(2,2) = mu(i)/3.0*(-eps*eps/l2 + 2.0*l2 - l3*l3/l2)/pow23 + mu(i)/9.0*(-4.0*eps/l2 - 4.0/eps*l2 + 2.0/eps/l2*l3*l3)/pow23*(l1-eps) + K*(Ja(i)-1.0)*l1*l3;
        Pd(3,3) = mu(i)/3.0*(-eps*eps/l3 - l2*l2/l3 + 2.0*l3)/pow23 + mu(i)/9.0*(-4.0*eps/l3 + 2.0/eps*l2*l2/l3 - 4.0/eps*l3)/pow23*(l1-eps) + K*(Ja(i)-1.0)*l1*l2;
%         W(i) = 0.5*mu(i)*((eps*eps + l2*l2 + l3*l3)/pow23 - 3.0) + mu(i)/3.0*(2.0*eps - l2*l2/eps - l3*l3/eps)/pow23*(l1-eps) + 0.5*k*(l1-eps)*(l1-eps) + 0.5*K*((J1(i)-1.0)*(J1(i)-1.0) + (J2(i)-1.0)*(J2(i)-1.0) + (J3(i)-1.0)*(J3(i)-1.0) + (J4(i)-1.0)*(J4(i)-1.0))/4.0;
        W(i)= rho;
        P{i} = U1*(Pd*U');
    end;
    
    % Apply forces to nodes
    Ft(:,El(i,1)) = Ft(:,El(i,1)) + P{i}*(N1{i} + N2{i} + N3{i})'/6.0;
    Ft(:,El(i,2)) = Ft(:,El(i,2)) + P{i}*(N1{i} + N3{i} + N4{i})'/6.0;
    Ft(:,El(i,3)) = Ft(:,El(i,3)) + P{i}*(N2{i} + N3{i} + N4{i})'/6.0;
    Ft(:,El(i,4)) = Ft(:,El(i,4)) + P{i}*(N1{i} + N2{i} + N4{i})'/6.0;	
    
    % Growth
    G{i} = eye(3)+(eye(3)-[Ns(i,1)*Ns(i,1),Ns(i,1)*Ns(i,2),Ns(i,1)*Ns(i,3);Ns(i,2)*Ns(i,1),Ns(i,2)*Ns(i,2),Ns(i,2)*Ns(i,3);Ns(i,3)*Ns(i,1),Ns(i,3)*Ns(i,2), Ns(i,3)*Ns(i,3)])*at*gm(i);
end;


% Midplane
mpy = -0.004;
for i = 1:nsn 
    pt = sn(i);
    if (V0(pt,2) < mpy-0.5*a) && (V(pt,2) > mpy)
        Ft(2,pt) = Ft(2,pt)-(mpy-V(pt,2))/hc*a*a*K;
    end;
    if (V0(pt,2) > mpy+0.5*a) && (V(pt,2) < mpy)   
        Ft(2,pt) = Ft(2,pt)-(mpy-V(pt,2))/hc*a*a*K;
    end;
end;

parfor i = 1:size(V,1)
    Ft(:,i) = Ft(:,i)-Vt(:,i)*gamma*Vn0(i);
    Vt(:,i) = Vt(:,i)+Ft(:,i)/(Vn0(i)*rho)*dt;
end;

t=t+dt;

% for i = 1:size(El,1)
% %     T{i} = (B{i}-eye(3)*trace(B{i})/3.0)*mu(i)/(J{i}*(J{i}^(2.0/3.0))) + eye(3)*K*(Ja(i)-1.0);
% %     P{i} = T{i}.*(inv((F{i})'))*J{i}; %for calculating the forces
% %     W{i} = 0.5*mu(i)*(trace(B{i})/J{i}^(2.0/3.0)-3.0)+0.5*K*((J1(i)-1.0)*(J1(i)-1.0)+(J2(i)-1.0)*(J2(i)-1.0)+(J3(i)-1.0)*(J3(i)-1.0)+(J4(i)-1.0)*(J4(i)-1.0))*0.25;
%     E(i) = (9*K*mu(i))/(3*K+mu(i));
%     v(i) = ((3*K) - (2*mu(i)))/((6*K) + (2*mu(i)));
% end;

% W(isnan(W) | isinf(W))=1;

% W = cell2mat(W);
% W=ones(size(W));
