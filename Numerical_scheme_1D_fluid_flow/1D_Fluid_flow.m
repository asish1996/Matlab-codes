

%% Common entries in all schemes
clear;
%tic;

u_boundaries = [40,0];
h = 0.04;

dy = 0.001;
num_nodes = h/dy + 1;
dt = 0.002;
%dt = 0.0025; %uncomment this to check Q7
v = 0.000217; % m^2/s
u_initial= zeros(num_nodes,1);
u_initial(1,1) = u_boundaries(1);  % 40m/s
u_initial(num_nodes,1) = u_boundaries(2);
u_final = zeros(num_nodes,1);
u_final (1,1) = u_boundaries(1);
u_final(num_nodes,1) = u_boundaries(2);
coeff = (v*dt)/(dy*dy);

u_intermediate = u_initial;

time = [0,0.18,0.36,0.54,0.72,0.9,1.08]; %plots required at different time steps
%time = [1.08]; %uncomment this to check validity of Q8
y = zeros(num_nodes,1);

for j = 2:num_nodes;
    y(j,1) = y(j-1,1) + dy;
end

%% FTCS Explicit scheme

u_explicit = [];

for t = 1:length(time);
    num_timesteps = time(t)/dt;
    
    for n = 1:num_timesteps;
        u_intermediate = u_final;
        
        for i = 2:num_nodes-1;
            u_final(i,1) = u_intermediate(i,1) + coeff * (u_intermediate(i+1,1) - 2*u_intermediate(i,1) + u_intermediate(i-1,1));
        end
    end
    u_explicit = [u_explicit,u_final]; %This array stores u at different nodes at different times as taken in time array
    %plot (u_final,y);
    %hold on;
    u_intermediate = u_initial;
    u_final = u_initial;
end

%toc;

%% FTCS Implicit scheme

a = -coeff;
b = 1+2*coeff;
c = -coeff;
u_implicit_intermediate = zeros(num_nodes-2,1);
u_implicit_final = zeros(num_nodes-2,1);
u_implicit_bc = zeros(num_nodes-2,1);
u_implicit_bc(1,1) = (-a)*u_boundaries(1);
u_implicit_bc(num_nodes-2,1) = (-c)*u_boundaries(2); 

u_implicit = [];
A_implicit = full(gallery('tridiag',num_nodes-2,a,b,c));

for t_implicit = 1:length(time);
    num_timesteps_implicit = time(t_implicit)/dt;
    
    for n_time = 1:num_timesteps_implicit;
        u_implicit_intermediate = u_implicit_final;        
        u_implicit_final = inv(A_implicit)*(u_implicit_intermediate + u_implicit_bc);
    end
    u_implicit_plot = [u_boundaries(1);u_implicit_final;u_boundaries(2)];
    u_implicit = [u_implicit,u_implicit_final]; %This array stores u at different nodes at different times as taken in time array
    %plot (u_implicit_plot,y);
    %hold on;
    u_implicit_intermediate = zeros(num_nodes-2,1);
    u_implicit_final = zeros(num_nodes-2,1);
end

first_row = ones(1,length(time)) * u_boundaries(1);
last_row = ones(1,length(time)) * u_boundaries(2);
u_implicit = [first_row;u_implicit;last_row];

%toc;

%% Crank Nicolson

A_cn = full(gallery('tridiag',num_nodes-2,-coeff,1+(2*coeff),-coeff));
u_cn_bc = zeros(num_nodes-2,1);
u_cn_bc(1,1) = coeff*u_boundaries(1);
u_cn_bc(num_nodes-2,1) = coeff*u_boundaries(2);

u_cn_intermediate = u_initial;
u_cn_final = u_initial;
u_cn = [];

for t_cn = 1:length(time);
    num_timesteps_cn = time(t_cn)/dt;
    
    for n_cn = 1:num_timesteps_cn;
        
        for beta = 2:num_nodes-1;
            u_cn_intermediate(beta) = coeff*u_cn_final(beta-1,1) + (1-2*coeff)*u_cn_final(beta,1) + coeff*u_cn_final(beta+1,1);
        end
        
        u_cn_final(2:num_nodes-1,1) = inv(A_cn)*(u_cn_intermediate(2:num_nodes-1,1) + u_cn_bc);
        %u_cn_final(2:num_nodes-1,1) = TDMAsolver(A_cn,(u_cn_intermediate(2:num_nodes-1,1) + u_cn_bc));
        %%uncomment the above step to verify TDMA approach
    end
    
    u_cn = [u_cn,u_cn_final];
    %plot (u_cn_final,y)
    %hold on
    
    u_cn_final = u_initial;
    u_cn_intermediate = u_initial;
    
end

%toc;

%% All methods at t=1.08                

%plot (u_explicit(:,length(time)),y)          %%Uncomment these to check Q6
%hold on

%plot (u_implicit(:,length(time)),y)
%hold on

%plot (u_cn(:,length(time)),y)


%% Computation time histogram
%%uncomment this to get bargraph of computation time

% methods = categorical({'FTCS Explicit', 'FTCS Implicit' , 'CN without TDMA' , 'CN with TDMA' });
% Elapsed_time = [0.00399;0.017611;0.015098;0.004511];  %noted separately for every method calculation by commenting out others
% bar(methods,Elapsed_time)


%% TDMA

function X=TDMAsolver(A,b)

m=length(b);                 % m is the number of rows
X=zeros(m,1);
A(1,2)= A(1,2)  ./ A(1,1);    % Division by zero risk.
b(1)=  b(1)    ./ A(1,1);    % Division by zero would imply a singular matrix

for i=2:m-1
    temp=  A(i,i) - A(i,i-1) .* A(i-1,i);
    A(i,i+1)=  A(i,i+1)  ./ temp;
    b(i)= ( b(i) - A(i,i-1) .* b(i-1) )  ./ temp;
end 

i=m;
X(m)=(b(i) - A(i,i-1) .* b(i-1))  ./ (A(i,i) - A(i,i-1) .* A(i-1,i));

for i=m-1:-1:1
X(i)=  -A(i,i+1) .* X(i+1) + b(i);
end
end





