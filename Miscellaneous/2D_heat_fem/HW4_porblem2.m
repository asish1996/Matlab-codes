% Shape functions for Q8
clc;
clear all;

% Input:
% Thermal conductivity

K = [10 0 ; 0 10];
% K = [1 0 ; 0 10];
%K = [1 0 ; 0 0.1];
%K = [3 1 ; 1 4];  % Chose the positive definite K 

nodes_per_elem = 4;     % nodes per element
no_elem_x = 50;          % no. of elements in the x-direction
no_elem_y =50;           % no. of elements in the y-direction
L_x= 2;                 %  length of rectangular domain in x-direction
L_y = 2;                % length of rectangular domian in y-direction
h_x=L_x/no_elem_x;         % step-size x -direction
h_y = L_y/no_elem_y;       % step-size y-direction

x_coord = linspace(0,L_x,no_elem_x+1);  
y_coord = linspace(0,L_y,no_elem_y+1);

% Boundary conditions
T_left=100;
T_right=200;
T_horz = @(x) (T_right-T_left)/L_x * x + T_left;        %  for top and bottom edges

x1=0;x2=1;x3=1;x4=0;
y1=0;y2=0;y3=1;y4=1;
syms xi eta %x1 x2 x3 x4 y1 y2 y3 y4

% shape functions for Q4
N1 = (1-xi)*(1-eta)/4;
N2 = (1+xi)*(1-eta)/4;
N3 = (1+xi)*(1+eta)/4;
N4 = (1-xi)*(1+eta)/4;

N = [N1 N2 N3 N4];

N_xi = diff(N,xi);
N_eta = diff(N,eta);

% Connectivity matrix
connectivity = zeros(no_elem_x*no_elem_y,4);
% intializing
x_loc = zeros(nodes_per_elem,1);
y_loc = zeros(nodes_per_elem,1);
K_global = zeros((no_elem_x+1)*(no_elem_y+1));


Dirichlet_boundary_bottom = [1:no_elem_x+1];
Dirichlet_boundary_right = [no_elem_x+1:(no_elem_x+1):(no_elem_x+1)*(no_elem_y+1)];
Dirichlet_boundary_top = [(no_elem_x+1)*(no_elem_y)+1:1:(no_elem_x+1)*(no_elem_y+1)];
Dirichlet_boundary_left = [1:no_elem_x+1:(no_elem_x+1)*no_elem_y+1];

solution_temp = zeros((no_elem_x+1)*(no_elem_y+1),1);

for i=1:size(solution_temp)
    if ismember(i,Dirichlet_boundary_left)
        solution_temp(i)=     T_left;
    end

    if ismember(i,Dirichlet_boundary_right)
        solution_temp(i)=     T_right;
    end

    if ismember(i,Dirichlet_boundary_bottom)
        solution_temp(i)=     T_horz((i-1)*h_x);
    end

    if ismember(i,Dirichlet_boundary_top)
        solution_temp(i)=     T_horz((i-Dirichlet_boundary_top(1))*h_x);
    end

end



for i=1:no_elem_x*no_elem_y
    r= ceil(i/no_elem_x);
    c= i-(r-1)*no_elem_x;
    ele_node_1 = (no_elem_x+1)* (r-1) + c;
    ele_node_2 = ele_node_1 +1;
    ele_node_3 = ele_node_2 + no_elem_x+1;
    ele_node_4 = ele_node_1 + no_elem_x+1;
    
    connectivity(i,1) = ele_node_1;
    connectivity(i,2) = ele_node_2;
    connectivity(i,3) = ele_node_3;
    connectivity(i,4) = ele_node_4;
    
    x_loc = [(c-1)*h_x  c*h_x c*h_x  (c-1)*h_x ];
   
    
    y_loc = [(r-1)*h_y  (r-1)*h_y r*h_y r*h_y];
    
    Jacobian = [N_xi; N_eta]* [x_loc' y_loc'];
Jaco_inv = inv(Jacobian);

R2=Jaco_inv;
R3 = [N_xi; N_eta];

B = R2*R3;


k_ele = int(int(B'*K*B*det(Jacobian),xi,-1,1),eta,-1,1);
    
    for m=1:nodes_per_elem
    for n=1:nodes_per_elem
        K_global(connectivity(i,m),connectivity(i,n)) = K_global(connectivity(i,m),connectivity(i,n)) + k_ele(m,n);
    end
    end

end
K_red= K_global;
load_global = zeros((no_elem_x+1)*(no_elem_y+1),1);
for i=size(solution_temp,1):-1:1
    if (solution_temp(i) ~= 0 )
        load_global = load_global - solution_temp(i)*K_global(:,i);
        
        K_red(i,:)=[];
        K_red(:,i) = [];
        
    end
end

for i=size(solution_temp,1):-1:1
    if (solution_temp(i) ~= 0)
        load_global(i)= [];
    end
end


T=inv(K_red)*load_global;
j=1;
for i=1:size(solution_temp,1)
    
    if (solution_temp(i) == 0)
        solution_temp(i) = T(j);
        j=j+1;
    end
    
end

%% Plot for element temperature

N_xi_ele = 10;
N_eta_ele = 10;
T_plot= zeros(N_xi_ele,N_eta_ele);
flux =zeros(no_elem_x* no_elem_y,2);
for i=1:no_elem_x*no_elem_y
N_xi_ele = 10;
N_eta_ele = 10;
T_plot= zeros(N_xi_ele,N_eta_ele);
flux =zeros(no_elem_x* no_elem_y,2);
for i=1:no_elem_x*no_elem_y
    r= ceil(i/no_elem_x);
    c= i-(r-1)*no_elem_x;
    ele_node_1 = (no_elem_x+1)* (r-1) + c;
    ele_node_2 = ele_node_1 +1;
    ele_node_3 = ele_node_2 + no_elem_x+1;
    ele_node_4 = ele_node_1 + no_elem_x+1;
       
    connectivity(i,1) = ele_node_1;
    connectivity(i,2) = ele_node_2;
    connectivity(i,3) = ele_node_3;
    connectivity(i,4) = ele_node_4;
       
    x_loc = [(c-1)*h_x  c*h_x c*h_x  (c-1)*h_x ];
     
    y_loc = [(r-1)*h_y  (r-1)*h_y r*h_y r*h_y];
    
    Jacobian = [N_xi; N_eta]* [x_loc' y_loc'];
Jaco_inv = inv(Jacobian);

R2=Jaco_inv;
R3 = [N_xi; N_eta];

B = R2*R3;
Temp_element = [solution_temp(ele_node_1), solution_temp(ele_node_2),solution_temp(ele_node_3), solution_temp(ele_node_4)]';

flux (i,:) = -K*B*Temp_element;
xi_vec = linspace(-1,1,N_xi_ele);
eta_vec = linspace(-1,1,N_eta_ele);
for m=1:N_xi_ele
    for n=1:N_eta_ele
        T_plot(m,n) = Temp_element'*subs(subs(N,xi,xi_vec(m)),eta,eta_vec(n))'; 
    end
end



end

    r= ceil(i/no_elem_x);
    c= i-(r-1)*no_elem_x;
    ele_node_1 = (no_elem_x+1)* (r-1) + c;
    ele_node_2 = ele_node_1 +1;
    ele_node_3 = ele_node_2 + no_elem_x+1;
    ele_node_4 = ele_node_1 + no_elem_x+1;
       
    connectivity(i,1) = ele_node_1;
    connectivity(i,2) = ele_node_2;
    connectivity(i,3) = ele_node_3;
    connectivity(i,4) = ele_node_4;
       
    x_loc = [(c-1)*h_x  c*h_x c*h_x  (c-1)*h_x ];
     
    y_loc = [(r-1)*h_y  (r-1)*h_y r*h_y r*h_y];
    
    Jacobian = [N_xi; N_eta]* [x_loc' y_loc'];
Jaco_inv = inv(Jacobian);

R2=Jaco_inv;
R3 = [N_xi; N_eta];

B = R2*R3;
Temp_element = [solution_temp(ele_node_1), solution_temp(ele_node_2),solution_temp(ele_node_3), solution_temp(ele_node_4)]';

flux (i,:) = -K*B*Temp_element;
xi_vec = linspace(-1,1,N_xi_ele);
eta_vec = linspace(-1,1,N_eta_ele);
for m=1:N_xi_ele
    for n=1:N_eta_ele
        T_plot(m,n) = Temp_element'*subs(subs(N,xi,xi_vec(m)),eta,eta_vec(n))'; 
    end
end



end
