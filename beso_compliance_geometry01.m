%% cleaning

clc; clear; close all

%% coordinate matrix

nx = 81; % nodes number in x
ny = 51; % nodes number in y
Lx = 80/1000; % length in x [mm]
Ly = 50/1000; % length in y [mm]
Lz = 1/1000; % length in z [mm]
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

%% boundary conditions matrix

bc = zeros(ny,3);
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

%% material matrix
            %E    rho  nu
material = [100e9 7000 0.3];  % [Pa] [kg/m�] []

%% solver

nload = size(load,1); %load quantity
nnodes = size(coord,1); % number of nodes
alldof = 1:nnodes*2; % degrees of freedom
nel = size(inci,1); % number of elements
F = zeros(2*nnodes,1); %load matrix pre-location 
kg = zeros(2*nnodes); % global stiffness matrix pre-location

for i = 1:nload % load vector
    if load(i,2) == 1
        F(2*load(i,1)-1) = load(i,3);        
    else
        F(2*load(i,1)) = load(i,3);
    end    
end

for i = 1 : nel % global stiffness matrix 

    no1 = inci(i,2); % first node element 
    no2 = inci(i,3); % second node element
    no3 = inci(i,4); % third node element
    no4 = inci(i,5); % fourth node element

    E = material(inci(i,1),1); % young's module

    nu = material(inci(i,1),3); % poisson module

    k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ... 
   -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];

    ke = Lz*E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
                      k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
                      k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
                      k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
                      k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
                      k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
                      k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
                      k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];

    loc = [no1*2-1 no1*2 no2*2-1 no2*2 no3*2-1 no3*2 no4*2-1 no4*2]; % localization vector
    
    kg(loc,loc) = kg(loc,loc) + ke; % global stiffness matrix assemble
 
end

freedof = alldof;

ke_min = 1e-3*ke;

for k = 1 : size(bc,1) % degrees of freedom
    freedof(2*bc(k,1)-(2-bc(k,2))) = 0;
end

F_aux = sparse(F(logical(freedof),1)); % column & rows elimination

u = zeros(size(F_aux,1),1);

%% BESO parameters

Vf = 0.5;
Er = 1/100;
Ar_max = 5/100;
r_min = 3/1000;
tau = 0.01/100;

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

figure()
patch('Faces',inci(:,2:5),'Vertices',coord(:,2:3),'FaceColor','blue')
axis equal
hold on
% chr = int2str(coord(:,1));
% text(coord(:,2),coord(:,3),chr)

%% BESO compliance optimization

presence = ones(nel,1); % elements presence

iter = 1 ; % iteration counter

nel_1 = nel; % volume fraction * n of full elements each iteraction

nel_f = nel * Vf; % n of elements at the end of the optimization

nel_add_max = nel * Ar_max; % max n of elements added each iteraction

complaince = zeros(10,1); % complaince value

historical = zeros(nel,2); % history of the sensivity number vector

elements = ones(nel,2); % presence of elements each two iteractions

complaince_hist = zeros(100,1);

conv = 0;

kk = 0;

tic

while conv == 0
    
    if nel_1 > nel_f
        
        nel_0 = nel_1;
    
        nel_1 = ceil(nel_1*(1-Er));
        
    end
       
    kg_aux = sparse(kg(logical(freedof),logical(freedof))); % column & rows elimination
    
    u(logical(freedof),1) = kg_aux\F_aux; % displacement vector
    
    alpha_elemental = zeros(nel,1); % sensivity number vector

    alpha_n = zeros(nnodes,1); % sensitivity nodal number
    
    C = 0; % complaince value 

    for i = 1 : nel % alpha_elemental / alpha_nodal
    
    no1 = inci(i,2); % first node element 
    no2 = inci(i,3); % second node element
    no3 = inci(i,4); % third node element
    no4 = inci(i,5); % fourth node element
    
    nodes = [no1,no2,no3,no4]; % set of nodes

    loc = [no1*2-1 no1*2 no2*2-1 no2*2 no3*2-1 no3*2 no4*2-1 no4*2]; % localization vector

    u_aux = u(loc,1); % elementar displacement vector
    
    if presence(i) == 1
        
        alpha_elemental(i,1) = 0.5*u_aux'*ke*u_aux; % sensivity vector
        
        C = C + alpha_elemental(i,1);
        
    else
        alpha_elemental(i,1) = 0.5*u_aux'*ke_min*u_aux; % sensivity vector
        
        C = C + alpha_elemental(i,1);
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
    
    complaince_hist(iter,1) = C;
    
    if iter <= 10 
        
        complaince(iter,1) = C;
    
    else
        
        iter_2 = mod(iter,10);
        
        if iter_2 == 0
            iter_2 = 10;
        end
        
        complaince(iter_2,1) = C;
        
        error = abs(sum(complaince(6:10)) - sum(complaince(1:5)))/sum(complaince(6:10));
        
    end
    
    alpha_i = zeros(nel,1);

    for i = 1 : nel 

        w_rij = nonzeros(distance(i,:));

        alpha_j = alpha_n(logical(node_set_aux(i,:)));

        alpha_i(i,1) = w_rij' * alpha_j/sum(w_rij);

    end % alpha_i
    
    historical(:,1) = historical(:,2);
    
    historical(:,2) = alpha_i;
    
    alpha_mean(:,1) = 1:nel;

    alpha_mean(:,2) = mean(historical,2);
    
    alpha_aux_th = sortrows(alpha_mean,-2);
    
    th = round(alpha_aux_th(nel_1,2),10);
    
    remove_elements = alpha_mean( (presence == 1) & alpha_mean(:,2)<= th);
    
    nel_del = length(remove_elements);
    
    add_elements = alpha_mean((presence == 0) & alpha_mean(:,2) > th);
    
    nel_add = length(add_elements);
    
    if nel_add > nel_add_max
        
        nel_add = nel_add_max;
        
        add_elements = alpha_mean(alpha_mean(presence == 0, 2) > th,:);
        
        add_elements = sortrows(add_elements, -2);
        
        add_elements = add_elements(1:nel_add);
        
        nel_del = nel_add + sum(presence)*Er;
        
        remove_elements = alpha_mean(alpha_mean(presence == 1, 2) <= th,:);
        
        remove_elements = sortrows(remove_elements, 2);
        
        remove_elements = remove_elements(1:nel_del);
        
    end
    
    presence(remove_elements(:,1)) = 0;
    
    presence(add_elements(:,1)) = 1;
    
     for j = 1:nel_del
         
        i = remove_elements(j,1);
         
        no1 = inci(i,2); % first node element 
        no2 = inci(i,3); % second node element
        no3 = inci(i,4); % third node element
        no4 = inci(i,5); % fourth node element

        loc = [no1*2-1 no1*2 no2*2-1 no2*2 no3*2-1 no3*2 no4*2-1 no4*2]; % localization vector

        kg(loc,loc) = kg(loc,loc) - ke + ke_min; % global stiffness matrix assemble
        
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

        kg(loc,loc) = kg(loc,loc) - ke_min + ke; % global stiffness matrix assemble
        
    end
     
    inci_aux = inci(add_elements(:,1),:);
    patch('Faces',inci_aux(:,2:5),'Vertices',coord(:,2:3),'FaceColor','blue')
    pause(1e-6);
    
    if nel_1 <= nel_f && error < tau
            
            conv = 1;
            
    end
    
    if mod(iter, 15) == 0
        kk = kk + 1;
        filename = sprintf('%d_%d_%d', ex, ey, kk);
        saveas(gcf,filename,'png')
    end
    
    iter = iter + 1;
    

end

toc

iter = 1:iter;

complaince_hist = complaince_hist(iter);

complaince_hist = nonzeros(complaince_hist);

iter = iter(1:length(complaince_hist));

figure('Name','Mean Complaince','NumberTitle','off');
plot(iter,complaince_hist)
ylabel('Mean complaince (Nmm)')
xlabel('Iteration')
