% UNIVERSIDADE ESTADUAL DE CAMPINAS
% FACULDADE DE ENGENHARIA MECANICA
% METODOS DE OTIMIZACAO TOPOLOGICA EVOLUCIONARIA - IM437 J
% ATIVIDADE 3
%
% DIPL. -ENG RENAN MIRANDA PORTELA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% cleaning

clc; clear; close all

%% coordinate matrix

nx = 65; % nodes number in x
ny = 41; % nodes number in y
Lx = 0.16; % length in x
Ly = 0.10; % length in y
ex = nx - 1; % element number in x
ey = ny - 1; % element number in y
a = Lx/nx; % element length
b = Ly/ny; % element hight

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

% figure()
% patch('Faces',inci(:,2:5),'Vertices',coord(:,2:3),'FaceColor','none')
% axis equal
% title('Otimização Topológica Evolucionária')
% ylabel('y')
% xlabel('x')
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

%% load matrix

if mod(ey,2) == 0
    k = nx*(ny+1)/2;
else
    k = nx*(ny)/2;
end

load = [k 2 -1000];

% figure()
% scatter(coord(:,2),coord(:,3),'x')
% axis equal
% chr = int2str(coord(:,1));
% text(coord(:,2),coord(:,3),chr)
% hold on
% scatter(coord(k,2),coord(k,3),'o','red')
% scatter(coord(bc(:,1),2),coord(bc(:,1),3),'o','green')

%% material matrix
            %E          nu
material = [210e9 7860 0.3;  %steel
            70e9 2700 0.27]; %aluminium
        
%% solver
        
nel = size(inci,1);    % element number
nnodes = size(coord,1);% nodes number

alldof = 1:nnodes*2; % degrees of freedom
kg = zeros(2*nnodes); % global stiffness matrix pre-location
F = zeros(2*nnodes,1); %load matrix pre-location 
nload = size(load,1); %load quantity

for i = 1:nload
    if load(i,2) == 1
        F(2*load(i,1)-1) = load(i,3);        
    else
        F(2*load(i,1)) = load(i,3);
    end    
end

%% global stiffness matrix assemble

for i = 1:nel 
    
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
            
        %Pontos de gaus
        np=2;
        pxi=[-1/sqrt(3),1/sqrt(3)];
        peta=[-1/sqrt(3),1/sqrt(3)];
        wi=1;
        wj=1;
        
        Nfun=4; %Nï¿½mero de nï¿½s do elemento
        ke=zeros(2*Nfun); %Inicializaï¿½ï¿½o da matriz elementar
        Bs=zeros(3,2*Nfun); %Inicializaï¿½ï¿½o da matriz de deformaï¿½ï¿½o
        
        for j=1:np
            xi=pxi(j);
            for k=1:np
                eta=peta(k);
                N=[1/4*(1-xi)*(1-eta), 1/4*(1+xi)*(1-eta), 1/4*(1+xi)*(1+eta), 1/4*(1-xi)*(1+eta)];
                dNxi=[eta/4 - 1/4,1/4 - eta/4,eta/4 + 1/4,- eta/4 - 1/4];
                dNeta=[xi/4 - 1/4,- xi/4 - 1/4,xi/4 + 1/4,1/4 - xi/4];
                %Jacobiano
                J=[dNxi;dNeta]*posxy;
                %Determinante do Jacobiano
                detJ=det(J);
                iJ=J^-1;
                %Matriz de deformaï¿½ï¿½o Bs-sï¿½lida
                p = 1:2:8; %vetor de 1 ate 8(
                Bs(1,p)   = iJ(1,1)*dNxi + iJ(1,2)*dNeta;
                Bs(2,p+1) = iJ(2,1)*dNxi + iJ(2,2)*dNeta;
                Bs(3,p)   = Bs(2,p+1);
                Bs(3,p+1) = Bs(1,p);
                %Matriz do elemento
                ke = ke + 0.001*Bs'*D*Bs*det(J)*wi*wj;
            end
        end
                              
        loc = [no1*2-1 no1*2 no2*2-1 no2*2 no3*2-1 no3*2 no4*2-1 no4*2]; % localization vector

        kg(loc,loc) = kg(loc,loc) + ke; % global stiffness matrix assemble
end

freedof = alldof;

for k = 1 : size(bc,1)
    freedof(2*bc(k,1)-(2-bc(k,2))) = 0;
end

F_aux = sparse(F(logical(freedof),1)); % column & rows elimination

u = zeros(size(F_aux,1),1);

figure()
patch('Faces',inci(:,2:5),'Vertices',coord(:,2:3),'FaceColor','blue')
axis equal
hold on

%% eso compliance optimization

tic

count = 0; % iteration count

err = 0.01; % element removal rate

nel_del = 0.01*nel; % quantity of element to be removed

nel_del = round(nel_del);

if mod(nel_del,2) ~= 0
    nel_del = nel_del - 1;
end

x=ones(nel,1); % elements presence

alpha_max = zeros(1000,1); % maximum sensitivity

u_max = zeros(1000,1); % maximum displacement

kk = 0;

nel_aux = zeros(1000,1);

while 1
    count = count + 1; % iteration count
    
%     nel_del = 0.002*sum(x); % quantity of element to be removed

    nel_del = 2;%round(nel_del)
    
    nel_aux(count,1) = sum(x);
    
%     if mod(nel_del,2) ~= 0
%         nel_del = nel_del + 1
%     end
    
    kg_aux = sparse(kg(logical(freedof),logical(freedof))); % column & rows elimination
    
    u(logical(freedof),1) = kg_aux\F_aux; % displacement vector
    
    u_max(count,1) = max(u);
    
    if max(u) >= 10e-5
        break
    end
    
    alpha = zeros(nel,2); % sensivity vector
    
    for i = 1 : nel
        
        no1 = inci(i,2); % first node element 
        no2 = inci(i,3); % second node element
        no3 = inci(i,4); % third node element
        no4 = inci(i,5); % fourth node element
        
        loc = [no1*2-1 no1*2 no2*2-1 no2*2 no3*2-1 no3*2 no4*2-1 no4*2]; % localization vector
               
        E = material(inci(i,1),1); % young's module

        nu = material(inci(i,1),3); % poisson module

        k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ... 
       -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];

        ke = 0.001*E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
                          k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
                          k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
                          k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
                          k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
                          k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
                          k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
                          k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
        
        u_aux = u(loc,1); % elementar displacement vector
        
        alpha(i,1) = i; % sensivity vector indice
        
        alpha(i,2) = 0.5*u_aux'*ke*u_aux; % sensivity vector
        
    end
    
    alpha_max(count,1) = max(alpha((x == 1),2));
    
    alpha = alpha((x == 1),:);
    
    alpha_aux = sortrows(alpha,2);
    
    remove_element_aux = alpha_aux(1:nel_del);
    
    x(remove_element_aux) = 0;
     
     for j = 1:nel_del
        i = remove_element_aux(j);
         
        no1 = inci(i,2); % first node element 
        no2 = inci(i,3); % second node element
        no3 = inci(i,4); % third node element
        no4 = inci(i,5); % fourth node element

        E = material(inci(i,1),1); % young's module

        nu = material(inci(i,1),3); % poisson module

        k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ... 
       -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];

        ke = 0.001*E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
                          k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
                          k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
                          k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
                          k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
                          k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
                          k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
                          k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];

        loc = [no1*2-1 no1*2 no2*2-1 no2*2 no3*2-1 no3*2 no4*2-1 no4*2]; % localization vector

        kg(loc,loc) = kg(loc,loc) - ke; % global stiffness matrix assemble
        
        ke = 1e-3*ke;
        
        kg(loc,loc) = kg(loc,loc) + ke; % global stiffness matrix assemble
     end
     
     inci_aux = inci(remove_element_aux,:);
     
     patch('Faces',inci_aux(:,2:5),'Vertices',coord(:,2:3),'FaceColor','white')
     pause(1e-6);
     
%     if mod(count, 3) == 0
%         kk = kk + 1;
%         filename = sprintf('%d_%d_%d', ex, ey, kk);
%         saveas(gcf,filename,'png')
%     end
end

toc

count = 1:count;

u_max = u_max(count,1);

figure()
plot(count,u_max,'k')
title('Deslocamento vs. iterações')
xlabel('Iterações')
ylabel('Deslocamento')

alpha_max = alpha_max(count,1);

figure()
plot(count,alpha_max,'k')
xlabel('Iterações')
ylabel('\alpha')

nel_aux = nel_aux(count,1);
figure()
plot(count,nel_aux,'k')
xlabel('Iterações')
ylabel('Número de elementos')