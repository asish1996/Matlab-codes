%%MATLAB CODE FOR FEA 2D TRUSS ELEMENT.

% 2D truss program                                                            

clc
clear all % clear memory
tic       % starts a stopwatch timer
% elements nodes
ele_nod=[1 2;1 3;2 3;2 4;1 4;3 4;3 6;4 5;4 6;3 5;5 6];
%number of elements
num_ele=size(ele_nod,1);
%number of nodes
num_nod=6;
% nodes coordinates
nod_coor=[0 0;0 3;3 0;3 3;6 0;6 3];
% elements degree of freedom (DOF) 
ele_dof=[1 2 3 4;1 2 5 6;3 4 5 6;...
3 4 7 8;1 2 7 8;5 6 7 8;5 6 11 12;...
7 8 9 10;7 8 11 12;5 6 9 10;9 10 11 12];
% A, E, L are cross sectional area, Young's modulus, length of elements,respectively.
A = ones(1,11);
A = 3*10^(-4) * A; %0.08333267
E = ones(1,11);
E = 70*10^9 * E; %36*10^7
% initial zero matrix for all matrices
displacement=zeros(2*num_nod,1);
force=zeros(2*num_nod,1);
stiffness=zeros(2*num_nod);
%applied loads at DOFs
force(4)=-1000000.0;
force(8)=-1000000.0;
force(11)=1000000.0;
force(12)=-1000000.0;
%Boundary conditions
displacement (1,1)=0.0;
displacement (2,1)=0.0;
displacement (10,1)=0.0;
% computation of the system stiffness matrix
for e=1:num_ele
 L(e)=sqrt((nod_coor(ele_nod(e,2),1)-nod_coor(ele_nod(e,1),1))^2+...
      (nod_coor(ele_nod(e,2),2)-nod_coor(ele_nod(e,1),2))^2);
 C=(nod_coor(ele_nod(e,2),1)-nod_coor(ele_nod(e,1),1))/L(e);
 S=(nod_coor(ele_nod(e,2),2)-nod_coor(ele_nod(e,1),2))/L(e);
 k=(A(e)*E(e)/L(e)*[C*C C*S -C*C -C*S;C*S S*S -C*S -S*S;...
   -C*C -C*S C*C C*S; -C*S -S*S C*S S*S]);
   
% extract the rows of ele_dof (for each element e)
ele_dof_vec=ele_dof(e,:);
    for i=1:4
        for j=1:4
  stiffness(ele_dof_vec(1,i),ele_dof_vec(1,j))=...
  stiffness(ele_dof_vec(1,i),ele_dof_vec(1,j))+k(i,j);
        end
    end
end
% known force array
known_f_a=[3;4;5;6;7;8;9;11;12];
for i=1:size(known_f_a,1)
   dis_new(i,1)=displacement(known_f_a(i,1),1);
   force_new(i,1)=force(known_f_a(i,1),1);
end
for i=1:size(known_f_a,1)
    for j=1:size(known_f_a,1)
    stiff_new(i,j)=stiffness(known_f_a(i,1),known_f_a(j,1));
end
end
% solving the partitioned matrix 
dis_new=stiff_new\force_new;
for i=1:size(known_f_a,1)
  displacement(known_f_a(i,1),1)=dis_new(i,1);
end
% undeformed truss plot
for e=1:num_ele
    x=[nod_coor(ele_nod(e,1),1) nod_coor(ele_nod(e,2),1)];
    y=[nod_coor(ele_nod(e,1),2) nod_coor(ele_nod(e,2),2)];
plot(x,y,'b')
hold on
end
% known dicplacement array
known_dis_a=[1;2;10];
for i=1:size(known_dis_a,1)
    force(known_dis_a(i,1),1)=stiffness(known_dis_a(i,1),:)...
    *displacement;
end
% stress in elements
for e=1:num_ele
 L(e)=sqrt((nod_coor(ele_nod(e,2),1)-nod_coor(ele_nod(e,1),1))^2+...
      (nod_coor(ele_nod(e,2),2)-nod_coor(ele_nod(e,1),2))^2);
 C=(nod_coor(ele_nod(e,2),1)-nod_coor(ele_nod(e,1),1))/L(e);
 S=(nod_coor(ele_nod(e,2),2)-nod_coor(ele_nod(e,1),2))/L(e);
 stress(e)=(E(e)/L(e))*[-C -S C S]*displacement((ele_dof(e,:))');
end
disp('stiffness')
stiffness
disp('displacement')
displacement
disp('force')
force
disp('stress')
stress'
k=0;
for i=1:6
    for j=1:2
    k=k+1;
    nod_coor_def(i,j)=nod_coor(i,j)+displacement(k,1);
    end
end
%deformed truss plot
for e=1:num_ele
    x=[nod_coor_def(ele_nod(e,1),1) nod_coor_def(ele_nod(e,2),1)];
    y=[nod_coor_def(ele_nod(e,1),2) nod_coor_def(ele_nod(e,2),2)];
plot(x,y,'r')
hold on
end
toc