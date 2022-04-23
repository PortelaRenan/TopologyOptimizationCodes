% @university: UNIVERSIDADE ESTADUAL DE CAMPINAS
% @school: FACULDADE DE ENGENHARIA MECANICA
% @module: METODOS DE OTIMIZACAO TOPOLOGICA EVOLUCIONARIA - IM437 J
% @activity: ATIVIDADE 6
%
% @author: DIPL. -ENG RENAN MIRANDA PORTELA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% cleaning

clc; clear; close all

%% coordinate matrix

nx = 321; % nodes number in x
ny = 41; % nodes number in y
Lx = 8; % length in x [m]
Ly = 1; % length in y [m]
Lz = 1/10; % length in z [m]
ex = nx - 1; % element number in x
ey = ny - 1; % element number in y
a = Lx/ex; 
b = Ly/ey;

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

%% boundary conditions matrix

bc(1,:) = [nx*(ny+1)/2-nx+1, 1, 0];
bc(2,:) = [nx*(ny+1)/2-nx+1, 2, 0];
bc(3,:) = [nx*(ny+1)/2, 1, 0];
bc(4,:) = [nx*(ny+1)/2, 2, 0];

%% material matrix
            %E    rho  nu
material = [10e6 1 0.3];  % [Pa] [kg/m2] []

%% code variables

nel = size(inci,1); % number of elements in the geometry
nnodes = size(coord,1); % number of nodes in the geometry

%% elemental stiffness/mass matrix 

E = material(inci(1,1),1); % young's module 

rho = material(inci(1,1),2); % density

nu = material(inci(1,1),3); % poisson module

D = E*[1/(1-nu^2), nu/(1-nu^2), 0; 
       nu/(1-nu^2),  1/(1-nu^2),0;
        0,       0, 1/(2*(1+nu))];

line = inci(1,:);

posxy = coord(line(2:length(line)-1),2:3);    

% Gauss Points
np= 2;
pxi=[-1/sqrt(3),1/sqrt(3)];
peta=[-1/sqrt(3),1/sqrt(3)];
wi=1;
wj=1;

Nfun=4; % Number of nodes by element
ke=zeros(2*Nfun); % Elemental stiffness matrix iniciation
Bs=zeros(3,2*Nfun); % Deformation matrix iniciation

for j=1:np
    xi=pxi(j);
    for k=1:np
        eta=peta(k);
        N=[1/4*(1-xi)*(1-eta), 1/4*(1+xi)*(1-eta), 1/4*(1+xi)*(1+eta), 1/4*(1-xi)*(1+eta)];
        dNxi=[eta/4 - 1/4,1/4 - eta/4,eta/4 + 1/4,- eta/4 - 1/4];
        dNeta=[xi/4 - 1/4,- xi/4 - 1/4,xi/4 + 1/4,1/4 - xi/4];
        %   Jacobian
        
        J=[dNxi;dNeta]*posxy;
        
        % Jacobian determination
        
        detJ=det(J);
        iJ=J^-1;
        
        % Deformation matrix
        
        p = 1:2:8; %vetor de 1 ate 8(
        Bs(1,p)   = iJ(1,1)*dNxi + iJ(1,2)*dNeta;
        Bs(2,p+1) = iJ(2,1)*dNxi + iJ(2,2)*dNeta;
        Bs(3,p)   = Bs(2,p+1);
        Bs(3,p+1) = Bs(1,p);
        
        %Elemental matrix
        
        ke = ke + Bs'*D*Bs*det(J)*wi*wj;
    end
end

%   Gauss points

PG=[-(1/sqrt(3)),-(1/sqrt(3));(1/sqrt(3)),-(1/sqrt(3));(1/sqrt(3)),(1/sqrt(3));-(1/sqrt(3)),(1/sqrt(3))];

me = zeros(8);

for j=1:4
    
    N=[(1/4)*(1-PG(j,1))*(1-PG(j,2)),0,(1/4)*(1+PG(j,1))*(1-PG(j,2)),0,(1/4)*(1+PG(j,1))*(1+PG(j,2)),0,(1/4)*(1-PG(j,1))*(1+PG(j,2)),0
        0,(1/4)*(1-PG(j,1))*(1-PG(j,2)),0,(1/4)*(1+PG(j,1))*(1-PG(j,2)),0,(1/4)*(1+PG(j,1))*(1+PG(j,2)),0,(1/4)*(1-PG(j,1))*(1+PG(j,2))];

    me=me+N'*N;
end

me = me*rho*a*b*0.25;
              
%% global stiffness/mass assembling 

kg = zeros(2*nnodes); % global stiffness matrix pre-location
mg = zeros(2*nnodes); % global mass matrix pre-location
alldof = 1:nnodes*2; % all degrees of freedom

for i = 1 : nel % global stiffness matrix 

    no1 = inci(i,2); % first node element 
    no2 = inci(i,3); % second node element
    no3 = inci(i,4); % third node element
    no4 = inci(i,5); % fourth node element
    
    loc = [no1*2-1 no1*2 no2*2-1 no2*2 no3*2-1 no3*2 no4*2-1 no4*2]; % localization vector
    
    kg(loc,loc) = kg(loc,loc) + ke; % global stiffness matrix assemble
    
    mg(loc,loc) = mg(loc,loc) + me; % global stiffness matrix assemble
    
end

freedof = alldof;

for k = 1 : size(bc,1) % free degrees of freedom
    freedof(2*bc(k,1)-(2-bc(k,2))) = 0;
end

%% BESO parameters

Vf = 0.5;
Er = 2/100; % Evolution rate
Ar_max = 5/100;
r_min = 0.075; % maximal distance
tau = 0.01/100;
p = 3; % penalization factor
x_min = 1e-6;
error = 100;

me_min = me*(1-x_min);
ke_min = ke*(1-x_min);

%% Filter

node_set = zeros(size(inci,1), 100);
node_set_aux = zeros(size(inci,1), size(coord,1));
distance = zeros(size(inci,1), 100);

for i = 1 : size(inci,1)
    
    x_nodes = coord(inci(i,2:5),2); % nodal coordinates in x
    y_nodes = coord(inci(i,2:5),3); % nodal coordinates in y

    x_c = mean(x_nodes,1); % coordinate of the center of the element in x
    y_c = mean(y_nodes,1); % coordinate of the center of the element in y
    
    k = 1;
    
    for j = 1 : size(coord,1)
        
        if (coord(j,2)-x_c)^2 + (coord(j,3)-y_c)^2 <= r_min^2
            node_set(i,k) = coord(j,1);
            node_set_aux(i,coord(j,1)) = 1;
            k = k + 1;
        end
        
    end
    
    for j = 1 : length(nonzeros(node_set(i,:)))
        
        distance(i,node_set(i,j)) = r_min - sqrt((coord(node_set(i,j),2)-x_c)^2 + (coord(node_set(i,j),3)-y_c)^2);
        
    end

    
end

figure('Name','Geometry','NumberTitle','off');
patch('Faces',inci(:,2:5),'Vertices',coord(:,2:3),'FaceColor','blue')
axis equal
hold on

%% BESO frequency optimization

presence = ones(nel,1); % elements presence
mode = 1; % mode 
iter = 0; % iteration counter
nel_1 = nel; % volume fraction * n of full elements each iteraction
nel_f = nel * Vf; % n of elements at the end of the optimization
nel_add_max = 10;%nel * Ar_max; % max n of elements added each iteraction
% historical = zeros(nel,2); % history of the sensivity number vector
historical = zeros(nel,50); % history of the sensivity number vector
frequency = zeros(1000,4);
Xi = ones(nel,1);
mean_C=zeros(1,10);
kk = 0;

conv = 0;

while conv == 0
    
    iter = iter + 1;
    
    if nel_1 > nel_f
        
        nel_0 = nel_1;
        nel_1 = ceil(nel_1*(1-Er));
        
        if nel_1 < nel_f
            
            nel_1 = nel_f;
            
        end
        
    else
        
        nel_1 = nel_f;
        
    end
    
    formatSpec = 'Numero de elementos na geometria: %4.0f/%4.0f \n';
    fprintf(formatSpec, nel_1, nel_f);
    
    kg_aux = sparse(kg(logical(freedof),logical(freedof))); % column & rows elimination
    mg_aux = sparse(mg(logical(freedof),logical(freedof))); % column & rows elimination
    [eig_vector,eig_value] = eigs(kg_aux,mg_aux,4,'sm'); % eigen_values & eigen_vectors
    
    eig_value = sum(eig_value);
    eig_aux = [eig_value;eig_vector]';
    eig_aux = sortrows(eig_aux,1);
    frequency(iter,:) = sqrt(eig_aux(:,1)); % 4 first frequencies
    omega = frequency(iter,mode); % mode of study frequency
    
    formatSpec = 'Omega: %4.2f \n';
    fprintf(formatSpec, omega);
    
    eig_vector = eig_aux(mode,2:end)';
    eig_vector_aux = zeros(size(kg,1),1);    
    eig_vector_aux(logical(freedof),1) = eig_vector;
    
    alpha_elemental = zeros(nel,1); % sensivity number vector
    alpha_n = zeros(nnodes,1); % sensitivity nodal number
    
    for i = 1 : nel % alpha_elemental / alpha_nodal
    
    no1 = inci(i,2); % first node element 
    no2 = inci(i,3); % second node element
    no3 = inci(i,4); % third node element
    no4 = inci(i,5); % fourth node element
    
    nodes = [no1,no2,no3,no4]; % set of nodes
    loc = [no1*2-1 no1*2 no2*2-1 no2*2 no3*2-1 no3*2 no4*2-1 no4*2]; % localization vector
    u_aux = eig_vector_aux(loc); % displacement nodes 
    
    if presence(i) == 1

        alpha_elemental(i,1) = (1/(2*omega))*u_aux'*(ke-(omega^2)/p*me)*u_aux;   

    else

        alpha_elemental(i,1) = -(omega/(2*p))*u_aux'*me*u_aux;

    end
        
        for j = 1:4

            x = coord(nodes(j),2); % nodal coordinate in x
            y = coord(nodes(j),3); % nodal coordinate in y

                if x == 0 || x == Lx 
                    if y == 0 || y == Ly
                        M = 1;
                    else
                        M = 2;
                    end
                elseif y == 0 || y == Ly
                    M = 2;
                else 
                    M = 4;
                end

            w = 1 / M;

            alpha_n(nodes(j)) = alpha_n(nodes(j)) + w * alpha_elemental(i);

        end
        
    end
    
    alpha_i = zeros(nel,1);

    for i = 1 : nel 

        w_rij = nonzeros(distance(i,:));
        alpha_j = alpha_n(logical(node_set_aux(i,:)));
        alpha_i(i,1) = w_rij' * alpha_j/sum(w_rij);

    end % alpha_i
    
    historical(:,iter) = alpha_i;
    
    alpha_mean(:,1) = 1:nel;
    alpha_mean(:,2) = mean(historical(:,1:iter),2);
    
%     historical(:,2) = alpha_i;
%     
%     alpha_mean(:,1) = 1:nel;
%     alpha_mean(:,2) = mean(historical,2);
%     historical(:,1) = alpha_mean(:,2);
    
    minimum = min(alpha_mean(:,2));
    maximum = max(alpha_mean(:,2));
    
    while abs((maximum-minimum)/maximum) > 1e-5
        
        th = (maximum + minimum)/2;
        Xi = max(x_min,sign(alpha_mean(:,2) - th));
        
        if sum(Xi) - nel_1 > 0
            minimum = th;
        else
            maximum = th;
        end
        
    end

    remove_elements = alpha_mean( (presence == 1) & alpha_mean(:,2) <= th);
    nel_del = length(remove_elements);
    
    add_elements = alpha_mean((presence == 0) & alpha_mean(:,2) > th);
    nel_add = length(add_elements);
    
    presence(remove_elements(:,1)) = 0;
    presence(add_elements(:,1)) = 1;
    
    for j = 1:nel_del
         
        i = remove_elements(j,1);
         
        no1 = inci(i,2); % first node element 
        no2 = inci(i,3); % second node element
        no3 = inci(i,4); % third node element
        no4 = inci(i,5); % fourth node element

        loc = [no1*2-1 no1*2 no2*2-1 no2*2 no3*2-1 no3*2 no4*2-1 no4*2]; % localization vector

        kg(loc,loc) = kg(loc,loc) - ke_min; % global stiffness matrix assemble
        mg(loc,loc) = mg(loc,loc) - me_min; % global stiffness matrix assemble

        
    end
    
    inci_aux = inci(remove_elements(:,1),:);
    patch('Faces',inci_aux(:,2:5),'Vertices',coord(:,2:3),'FaceColor','white')
    pause(1e-6);
    
    for j = 1:nel_add
         
        i = add_elements(j,1);
         
        no1 = inci(i,2); % first node element 
        no2 = inci(i,3); % second node element
        no3 = inci(i,4); % third node element
        no4 = inci(i,5); % fourth node element

        loc = [no1*2-1 no1*2 no2*2-1 no2*2 no3*2-1 no3*2 no4*2-1 no4*2]; % localization vector

        kg(loc,loc) = kg(loc,loc) + ke_min; % global stiffness matrix assemble
        mg(loc,loc) = mg(loc,loc) + me_min; % global stiffness matrix assemble
        
    end
     
    inci_aux = inci(add_elements(:,1),:);
    patch('Faces',inci_aux(:,2:5),'Vertices',coord(:,2:3),'FaceColor','blue')
    pause(1e-6);
    
    if iter < 10
        mean_C(1,iter)=omega;
    else
        mean_C(:,1)=[];
        mean_C(1,10)=omega;
        error=abs((sum(mean_C(1,6:10))-sum(mean_C(1,1:5)))/(sum(mean_C(1,6:10))));
    end
    
%     if mod(iter, 5) == 0
%         kk = kk + 1;
%         filename = sprintf('atividade_6_%d_%d_%d', ex, ey, kk);
%         saveas(gcf,filename,'fig')
%         saveas(gcf,filename,'png')
%     end
    
    if error < tau && nel_1 <= nel_f
        
        conv = 1;
        
    end
    
    formatSpec = 'Erro: %4.4f \n';
    fprintf(formatSpec, error);
    
end

iter = 1:iter;
frequency = frequency(iter,:);

frequency_1 = frequency(:,1);
frequency_2 = frequency(:,2);
frequency_3 = frequency(:,3);
frequency_4 = frequency(:,4);

figure('Name','Frequencies','NumberTitle','off');
hold on
grid on
grid minor

plot(iter,frequency_1,'o-')
plot(iter,frequency_2,'*-')
plot(iter,frequency_3,'^-')
plot(iter,frequency_4,'x')

% filename = sprintf('atividade_6_frequencies');
% saveas(gcf,filename,'fig')
% saveas(gcf,filename,'png')