% UNIVERSIDADE ESTADUAL DE CAMPINAS
% FACULDADE DE ENGENHARIA MECANICA
% METODOS DE OTIMIZACAO TOPOLOGICA EVOLUCIONARIA - IM437 J
% ATIVIDADE 7
%
% DIPL. -ENG RENAN MIRANDA PORTELA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% cleaning

clc; clear; close all

%% coordinate matrix

nx = 21; % nodes number in x
Lx = 20; % length in x [m]
ex = nx - 1; % element number in x
L = Lx/ex; 

coordx = linspace(0,Lx,nx); % coordinate in x
num = 1:nx; % number of nodes
coord = [num(:),coordx(:)]; % coordinate matrix

%% code variables

m_disc = 10;
I_disc = 100;
node_disc = (ex/2)+1;

%% incidence matrix

inci = zeros(ex,4); % incidence matrix pre-location

for i = 1 : ex
    
    inci(i,:) = [i coord(i,1) coord(i+1,1) 6];
    
end

%% boundary conditions matrix

bc(1,:) = [ coord(1,1), 1, 0];
bc(2,:) = [ coord(1,1), 2, 0];
bc(3,:) = [ coord(nx,1), 1, 0];
bc(4,:) = [ coord(nx,1), 2, 0];

%% material matrix
            %E    rho  nu
material = [2e11 7800 0.3];  % [Pa] [kg/m3] []

%% geometry

area = linspace(0.5, 1.5, 11);

%% solver

nnodes = size(coord,1); % number of nodes
alldof = 1:nnodes*2; % degrees of freedom
nel = size(inci,1); % number of elements
kg = zeros(2*nnodes); % global stiffness matrix pre-location
mg = zeros(2*nnodes); % global mass matrix pre-location

figure('Name','Geometry','NumberTitle','off');
xlim([0 L])
ylim([-5 5])
hold on
d = area(inci(:,4));
step=0;
for i=1:nel
    rectangle('Position',[0+step -d(i) L/nel 2*d(i)],'FaceColor','b')
    step=step+L/nel;
end

for i = 1 : nel % global stiffness matrix 

    no1 = inci(i,2); % first node element 
    no2 = inci(i,3); % second node element

    E = material(1,1); % young's module
    
    rho = material(1,2); % density
    
    A = area(inci(i,4))*0.1;
    
    I = A/(4*pi);

    ke = E*I/(L^3)*[12 6*L -12 6*L;
                    6*L 4*L^2 -6*L 2*L^2;
                    -12 -6*L 12 -6*L;
                    6*L 2*L^2 -6*L 4*L^2];
                      
    me = rho*A*L/420*[156 22*L 54 -13*L;
                        22*L 4*L^2 13*L -3*L^2;
                        54 13*L 156 -22*L;
                        -13*L -3*L^2 -22*L 4*L^2];

    loc = [no1*2-1 no1*2 no2*2-1 no2*2]; % localization vector
    
    kg(loc,loc) = kg(loc,loc) + ke; % global stiffness matrix assemble
    
    mg(loc,loc) = mg(loc,loc) + me; % global stiffness matrix assemble
 
end

freedof = alldof;

for k = 1 : size(bc,1) % degrees of freedom
    freedof(2*bc(k,1)-(2-bc(k,2))) = 0;
end

%% Adding mass to the mass matrix

loc_node_disc = [2*node_disc-1 2*node_disc];

mg(loc_node_disc(1),loc_node_disc(1)) = mg(loc_node_disc(1),loc_node_disc(1)) + m_disc; % global stiffness matrix assemble
mg(loc_node_disc(2),loc_node_disc(2)) = mg(loc_node_disc(2),loc_node_disc(2)) + I_disc; % global stiffness matrix assemble

%% BESO parameters

tau = 0.1/100;
mode = 1; % mode 
iter = 0 ; % iteration counter
frequency = zeros(100,4);
conv = 0;

while conv == 0 
    
    iter = iter + 1;

    kg_aux = (kg(logical(freedof),logical(freedof))); % column & rows elimination
    mg_aux = (mg(logical(freedof),logical(freedof))); % column & rows elimination
    
    [eig_vector,eig_value] = eig(kg_aux,mg_aux); %eig_vector = vibration mode
    eig_value = sum(eig_value);
    eig_aux = [eig_value;eig_vector]';
    omega = sqrt(eig_aux(1:4,1))/(2*pi);
    frequency(iter,:) = omega;
 
    eig_vector = eig_aux(mode,2:end)';
    eig_vector_aux = zeros(size(kg,1),1);   
    eig_vector_aux(logical(freedof),1) = eig_vector;
    
    alpha_elemental = zeros(nel,2); % sensivity number vector
    alpha_elemental(:,1) = 1:nel; 

    for i = 1 : nel
        
        no1 = inci(i,2); % first node element 
        no2 = inci(i,3); % second node element
        
        loc = [no1*2-1 no1*2 no2*2-1 no2*2]; % localization vector
        
        u_aux = eig_vector_aux(loc); 
        A = area(inci(i,4))*0.1;
        I = A/(4*pi);

        ke_i = E*I/(L^3)*[12 6*L -12 6*L;
                        6*L 4*L^2 -6*L 2*L^2;
                        -12 -6*L 12 -6*L;
                        6*L 2*L^2 -6*L 4*L^2];

        me_i = rho*A*L/420*[156 22*L 54 -13*L;
                            22*L 4*L^2 13*L -3*L^2;
                            54 13*L 156 -22*L;
                            -13*L -3*L^2 -22*L 4*L^2];
        
        alpha_elemental(i,2) = -u_aux'*(eig_value(mode)*me_i-ke_i)*u_aux;

    end
    
    alpha_aux = sortrows(alpha_elemental,-2);
    
    bigger = alpha_aux(1:2,:);
    
    for i = 1 : 2
        
        if inci(bigger(i,1),4) <= 10
        
            inci(bigger(i,1),4) = inci(bigger(i,1),4) + 1;

            no1 = inci(bigger(i,1),2); % first node element 
            no2 = inci(bigger(i,1),3); % second node element

            loc = [no1*2-1 no1*2 no2*2-1 no2*2]; % localization vector

            A = area(inci(bigger(i,1),4))*0.1;

            I = A/(4*pi);

            ke_i = E*I/(L^3)*[12 6*L -12 6*L;
                            6*L 4*L^2 -6*L 2*L^2;
                            -12 -6*L 12 -6*L;
                            6*L 2*L^2 -6*L 4*L^2];

            me_i = rho*A*L/420*[156 22*L 54 -13*L;
                                22*L 4*L^2 13*L -3*L^2;
                                54 13*L 156 -22*L;
                                -13*L -3*L^2 -22*L 4*L^2];

            kg(loc,loc) = kg(loc,loc) - ke + ke_i; % global stiffness matrix assemble

            mg(loc,loc) = mg(loc,loc) - me + me_i; % global stiffness matrix assemble
            
        end
        
    end
    
    if sum(area(inci(:,4))) > 20
    
        lower = alpha_aux(nel-1:nel,:);

        for i = 1 : 2

            if inci(lower(i,1),4) > 1

                inci(lower(i,1),4) = inci(lower(i,1),4) - 1;

                no1 = inci(lower(i,1),2); % first node element 
                no2 = inci(lower(i,1),3); % second node element

                loc = [no1*2-1 no1*2 no2*2-1 no2*2]; % localization vector

                A = area(inci(lower(i,1),4))*0.1;

                I = A/(4*pi);

                ke_i = E*I/(L^3)*[12 6*L -12 6*L;
                                6*L 4*L^2 -6*L 2*L^2;
                                -12 -6*L 12 -6*L;
                                6*L 2*L^2 -6*L 4*L^2];

                me_i = rho*A*L/420*[156 22*L 54 -13*L;
                                    22*L 4*L^2 13*L -3*L^2;
                                    54 13*L 156 -22*L;
                                    -13*L -3*L^2 -22*L 4*L^2];

                kg(loc,loc) = kg(loc,loc) - ke + ke_i; % global stiffness matrix assemble

                mg(loc,loc) = mg(loc,loc) - me + me_i; % global stiffness matrix assemble

            end

        end
        
    end
    
    if iter > 1
        
        error = (frequency(iter, mode) - frequency(iter-1, mode))/frequency(iter, mode);
        
        if error < tau 
            
            conv = 1;
            
        end
        
    end
    
end

% figure('Name','Geometry','NumberTitle','off');
% xlim([0 L])
% ylim([-5 5])
hold on
d = area(inci(:,4));
step=0;
color=gray(length(area));
for i=1:nel
    color1=color(inci(i,4),:);
    rectangle('Position',[0+step -d(i) L/nel 2*d(i)],'FaceColor',color1)
    step=step+L/nel;
end

hold off

iter = 1 : iter;
frequency = frequency(iter,:);

figure('Name','Frequencies','NumberTitle','off');
hold on
plot(frequency(:,1),'-o')
plot(frequency(:,2),'-^')
plot(frequency(:,3),'k-x')
plot(frequency(:,4),'g--')
xlabel('Iterations') 
ylabel('Frequency (Hz)')
legend({'1st Frequency','2nd Frequency','3rd Frequency','4th Frequency'})
title('Frequencies by iterations')
grid on
grid minor