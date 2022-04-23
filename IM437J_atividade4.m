% @university: UNIVERSIDADE ESTADUAL DE CAMPINAS
% @school: FACULDADE DE ENGENHARIA MECANICA
% @module: METODOS DE OTIMIZACAO TOPOLOGICA EVOLUCIONARIA - IM437 J
% @activity: ATIVIDADE 4
%
% @author: DIPL. -ENG RENAN MIRANDA PORTELA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% cleaning

clc; clear; close all

%% coordinate matrix

nx = 51; % nodes number in x
ny = 11; % nodes number in y
Lx = 5; % length in x [m]
Ly = 1; % length in y [m]
ex = nx - 1; % element number in x
ey = ny - 1; % element number in y

coordx = linspace(0,Lx,nx); % coordinate in x
coordy = linspace(0,Ly,ny); % coordinate in y
[X,Y] = meshgrid(coordx,coordy);
X = X';
Y = Y';
num = 1:nx*ny; % number of nodes
coord = [num(:),X(:),Y(:)]; % coordinate matrix

%% incidence matrix

inci = zeros(ex*ey,6); % incidence matrix pre-location
A = coord(:,1);
A = reshape(A,[nx,ny]);
k = 1;
for j = 1:ny-1
    for i = 1:nx-1
        if k > ex*ey/2
            mat = 1;
        else
            mat = 1;
        end
        inci(k,:)=[mat, A(i,j),A(i+1,j),A(i+1,j+1),A(i,j+1),k];
        k = k + 1;
    end
end

%% plot

figure()
patch('Faces',inci(:,2:5),'Vertices',coord(:,2:3),'FaceColor','blue')
axis equal
title('Otimização Topológica Evolucionária')
ylabel('y')
xlabel('x')
% chr = int2str(coord(:,1));
% text(coord(:,2),coord(:,3),chr)

%% boundary conditions matrix

bc = zeros(2*ny,3);
k = 1;

for i = 1:ny
    bc(2*i-1,:) = [k, 1,0];
    bc(2*i,:) = [k, 2,0];
    k = k + nx;
end

k = nx;

for i = ny+1:2*ny
    bc(2*i-1,:) = [k, 1,0];
    bc(2*i,:) = [k, 2,0];
    k = k + nx;
end

%% plot of nodes

% figure()
% scatter(coord(:,2),coord(:,3),'x')
% axis equal
% hold on
% scatter(coord(bc(:,1),2),coord(bc(:,1),3),'o','green')

%% material matrix
            %E    rho  nu
material = [200e9 7000 0.3;  %steel [Pa] [kg/m³] []
            70e9 2700 0.27]; %aluminium
        
%% solver
        
nel = size(inci,1);    % element number
nnodes = size(coord,1);% nodes number

alldof = 1:nnodes*2; % degrees of freedom
kg = zeros(2*nnodes); % global stiffness matrix pre-location
mg = zeros(2*nnodes); % global mass matrix pre-location

for i = 1:nel 
    
        no1 = inci(i,2); % first node element 
        no2 = inci(i,3); % second node element
        no3 = inci(i,4); % third node element
        no4 = inci(i,5); % fourth node element

        E = material(inci(i,1),1); % young's module
        
        rho = material(inci(i,1),2); % density

        nu = material(inci(i,1),3); % poisson module

        D = E*[1/(1-nu^2), nu/(1-nu^2), 0; 
               nu/(1-nu^2),  1/(1-nu^2),0;
                0,       0, 1/(2*(1+nu))];
        
        line = inci(1,:);
                
        posxy = coord(line(2:length(line)-1),2:3);    
        
        % Gauss points
        
        np=2;
        pxi=[-1/sqrt(3),1/sqrt(3)];
        peta=[-1/sqrt(3),1/sqrt(3)];
        wi=1;
        wj=1;
        
        Nfun=4; 
        ke=zeros(2*Nfun); 
        Bs=zeros(3,2*Nfun);
        
        teste = zeros(4);
        teste2 = zeros(4);
        teste3 = zeros(4);
        
        for j=1:np
            xi=pxi(j);
            for k=1:np
                eta=peta(k);
                N=[1/4*(1-xi)*(1-eta), 1/4*(1+xi)*(1-eta), 1/4*(1+xi)*(1+eta), 1/4*(1-xi)*(1+eta)];
                dNxi=[eta/4 - 1/4,1/4 - eta/4,eta/4 + 1/4,- eta/4 - 1/4];
                dNeta=[xi/4 - 1/4,- xi/4 - 1/4,xi/4 + 1/4,1/4 - xi/4];
                
                J=[dNxi;dNeta]*posxy;
                
                detJ=det(J);
                iJ=J^-1;
                
                p = 1:2:8; 
                B = iJ*[dNxi;dNeta];
                Bs(1,p)   = iJ(1,1)*dNxi + iJ(1,2)*dNeta;
                Bs(2,p+1) = iJ(2,1)*dNxi + iJ(2,2)*dNeta;
                Bs(3,p)   = Bs(2,p+1);
                Bs(3,p+1) = Bs(1,p);
                
                ke = ke + 0.1*Bs'*D*Bs*det(J)*wi*wj; % Elementar stiffness matrix
                teste = teste + B(1,:)'*B(1,:)*detJ;
                teste2 = teste2 + B(2,:)'*B(2,:)*detJ;
                teste3 = teste3 + N'*N*detJ;
            end
        end
        
        a = coord(no2,2) - coord(no1,2);
        
        b = coord(no3,3) - coord(no2,3);
        
        Area = a*b;
        
        me = 0.1*rho*Area/36*[4 0 2 0 1 0 2 0;
                              0 4 0 2 0 1 0 2;
                              2 0 4 0 2 0 1 0;
                              0 2 0 4 0 2 0 1;
                              1 0 2 0 4 0 2 0;
                              0 1 0 2 0 4 0 2;
                              2 0 1 0 2 0 4 0;
                              0 2 0 1 0 2 0 4];
                              
        loc = [no1*2-1 no1*2 no2*2-1 no2*2 no3*2-1 no3*2 no4*2-1 no4*2]; % localization vector

        kg(loc,loc) = kg(loc,loc) + ke; % global stiffness matrix assemble
        
        mg(loc,loc) = mg(loc,loc) + me; % global stiffness matrix assemble
end

freedof = alldof';

for k = 1 : size(bc,1)
    freedof(2*bc(k,1)-(2-bc(k,2))) = 0;
end

%% eso frequency optimization

count = 0; % iteration count

x=ones(nel,1); % elements presence

mode = 1; % vibration mode

alpha_max = zeros(1000,1);

frequency = zeros(1000,4);

nel_del = 4;

while 1
    
    count = count + 1; % iteration count
    
    kg_aux = sparse(kg(logical(freedof),logical(freedof))); % column & rows elimination
    
    mg_aux = sparse(mg(logical(freedof),logical(freedof))); % column & rows elimination
    
    [eig_vector,eig_value] = eigs(kg_aux,mg_aux,30,'SM');
    
    eig_value = sum(eig_value);
    
    eig_aux = [eig_value;eig_vector]';
    
    eig_aux = sortrows(eig_aux,1);
    
    omega = eig_aux(mode,1);
    
    eig_vector = eig_aux(mode,2:end)';
    
    eig_vector_aux = zeros(size(kg,1),1);
    
    eig_vector_aux(logical(freedof),1) = eig_vector;
    
    alpha = zeros(nel,2);
    
    frequency(count,1) = sqrt(omega)/(2*pi);
    
    frequency(count,2) = sqrt(eig_aux(2,1))/(2*pi);
    
    frequency(count,3) = sqrt(eig_aux(3,1))/(2*pi);
    
    frequency(count,4) = sqrt(eig_aux(4,1))/(2*pi);
    
    for i = 1:nel
        
        no1 = inci(i,2); % first node element 
        no2 = inci(i,3); % second node element
        no3 = inci(i,4); % third node element
        no4 = inci(i,5); % fourth node element
        
        element = [no1*2-1 no1*2 no2*2-1 no2*2 no3*2-1 no3*2 no4*2-1 no4*2];
        
        u = eig_vector_aux(element,1);      
        
        alpha(i,1) = i;
        
        alpha(i,2) = u'*(omega*me-ke)*u;               
        
    end
    
    if sum(alpha(:,2)) ~= 0
        print('Soma diferente de zero')
    end
    
    alpha_max(count,1) = max(alpha((x == 1),2));
    
    alpha = alpha((x == 1),:);
    
    alpha_aux = sortrows(alpha,-2);

    remove_element_aux = alpha_aux(1:nel_del,1);

    x(remove_element_aux) = 0;

    inci_aux = inci(remove_element_aux,:);

    patch('Faces',inci_aux(:,2:5),'Vertices',coord(:,2:3),'FaceColor','white')
    pause(1e-7)
    
     for j = 1:nel_del
        i = remove_element_aux(j);
         
        no1 = inci(i,2); % first node element 
        no2 = inci(i,3); % second node element
        no3 = inci(i,4); % third node element
        no4 = inci(i,5); % fourth node element

        E = material(inci(i,1),1); % young's module

        nu = material(inci(i,1),3); % poisson module

        D = E*[1/(1-nu^2), nu/(1-nu^2), 0; 
               nu/(1-nu^2),  1/(1-nu^2),0;
                0,       0, 1/(2*(1+nu))];
        
        line = inci(1,:);
                
        posxy = coord(line(2:length(line)-1),2:3);    
        
        % Gauss points
        
        np=2;
        pxi=[-1/sqrt(3),1/sqrt(3)];
        peta=[-1/sqrt(3),1/sqrt(3)];
        wi=1;
        wj=1;
        
        Nfun=4; 
        ke=zeros(2*Nfun); 
        Bs=zeros(3,2*Nfun);
        
        for j=1:np
            xi=pxi(j);
            for k=1:np
                eta=peta(k);
                N=[1/4*(1-xi)*(1-eta), 1/4*(1+xi)*(1-eta), 1/4*(1+xi)*(1+eta), 1/4*(1-xi)*(1+eta)];
                dNxi=[eta/4 - 1/4,1/4 - eta/4,eta/4 + 1/4,- eta/4 - 1/4];
                dNeta=[xi/4 - 1/4,- xi/4 - 1/4,xi/4 + 1/4,1/4 - xi/4];
                
                J=[dNxi;dNeta]*posxy;
                
                detJ=det(J);
                iJ=J^-1;
                
                p = 1:2:8; 
                Bs(1,p)   = iJ(1,1)*dNxi + iJ(1,2)*dNeta;
                Bs(2,p+1) = iJ(2,1)*dNxi + iJ(2,2)*dNeta;
                Bs(3,p)   = Bs(2,p+1);
                Bs(3,p+1) = Bs(1,p);
                
                ke = ke + 0.1*Bs'*D*Bs*det(J)*wi*wj; % Elementar stiffness matrix
            end
        end

        loc = [no1*2-1 no1*2 no2*2-1 no2*2 no3*2-1 no3*2 no4*2-1 no4*2]; % localization vector

        kg(loc,loc) = kg(loc,loc) - ke; % global stiffness matrix assemble
        
        ke = 1e-3*ke;
        
        kg(loc,loc) = kg(loc,loc) + ke; % global stiffness matrix assemble
        
        me = 0.1*rho*Area/36*[4 0 2 0 1 0 2 0;
                          0 4 0 2 0 1 0 2;
                          2 0 4 0 2 0 1 0;
                          0 2 0 4 0 2 0 1;
                          1 0 2 0 4 0 2 0;
                          0 1 0 2 0 4 0 2;
                          2 0 1 0 2 0 4 0;
                          0 2 0 1 0 2 0 4];
        
        mg(loc,loc) = mg(loc,loc) - me; % global stiffness matrix assemble
        
        me = 1e-3*me;
        
        mg(loc,loc) = mg(loc,loc) + me; % global stiffness matrix assemble
     end
    
    if count == 50
        break
    end
    
end

figure()
hold on
count = 1:count;
frequency_1 = frequency(count,1);
plot(count,frequency_1,'r');

frequency_2 = frequency(count,2);
plot(count,frequency_2,'b');

frequency_3 = frequency(count,3);
plot(count,frequency_3,'k');

frequency_4 = frequency(count,4);
plot(count,frequency_4,'g');

legend('\Omega_{1}','\Omega_{2}','\Omega_{3}','\Omega_{4}')
