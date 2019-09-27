load Water_Temp_Strain.txt

sensors = Water_Temp_Strain(:,:);
W = sensors(:,1); %Water level in meters
T = sensors(:,2); %Temperature in celsius
S = sensors(:,3); %Strain in micro starins
S = S %*10^-6;

%% a,b,c parts

nrows = [1:3650];
rowindex = [0;nrows'];
delT = 144;
time = rowindex.*delT; 
beta_mat = ones(3651);
beta_mat = beta_mat(:,1);

A = [time beta_mat T W];
c_m = (A'*A)^-1*A' * S;

q = c_m(1,1)*time + c_m(2,1)*beta_mat;
y = c_m(1,1)*time + c_m(2,1)*beta_mat + c_m(3,1)*T + c_m(4,1)*W;
q(1500,1)   % q1500
y(1500,1)   % strain at t1500

q3y = c_m(1,1)*(3650*3*144) + c_m(2,1) %actual strain after 3years
 %plot (time,q,'r')
%  hold on
%plot (time,y,'g')

%% d part

r = S - q; %residuals
%plot (time,r,'b')
d1 = mean(r)
d2 = std(r)

%% Chauvenet's criterion

 sy = std(r);
 my = mean(r);
 cs = [];
 cq = [];
 el = [];
 ctime = [];
 cW = [];
 cT = [];
for i = 1:3651;
    %range = [i];
        if r(i,1) < (my + 3*sy) && r(i,1) > (my - 3*sy);            
             
              ctime = [ctime;time(i,1)];
              cs = [cs;S(i,1)];
              cq = [cq;q(i,1)];
              cW = [cW;W(i,1)];
              cT = [cT;T(i,1)];
        else;
            el = [el;time(i,1)];
        end
        
%     end
end

rn1 = cs - cq;
sz2 = size(rn1);
cs2 = [];
cq2 = [];
el2 = [];
ct2 = [];
cW2 = [];
cT2 = [];

for i=1:sz2(1,1);
    if rn1(i,1) < (mean(rn1) + 3*std(rn1)) && rn1(i,1) > (mean(rn1) - 3*std(rn1));
        cs2 = [cs2;cs(i,1)];
        cq2 = [cq2;cq(i,1)];
        ct2 = [ct2;ctime(i,1)];
        cW2 = [cW2;cW(i,1)];
        cT2 = [cT2;cT(i,1)];
    else;
        el2 = [el2;ctime(i,1)];
    end
end

rn2 = cs2 - cq2;
sz3 = size(rn2);
cs3 = [];
cq3 = [];
el3 = [];
ct3 = [];
cW3 = [];
cT3 = [];

for i=1:sz3(1,1);
    if rn2(i,1) < (mean(rn2) + 3*std(rn2)) && rn2(i,1) > (mean(rn2) - 3*std(rn2));
        cs3 = [cs3;cs2(i,1)];
        cq3 = [cq3;cq2(i,1)];
        ct3 = [ct3;ct2(i,1)];
        cW3 = [cW3;cW2(i,1)];
        cT3 = [cT3;cT2(i,1)];
    else;
        el3 = [el3;ct2(i,1)];
    end
end

% We can see that cs2 = cs3 which indicates that solution has converged in
% 3 iterations


 %% With cleaned data

sz = size(cs3)
beta_mat2 = ones(sz(1,1),1);
a = [ct3, beta_mat2, cT3, cW3];
c_m2 = (a'*a)^-1*a' * cs3; %new w^

q2 = c_m2(1,1)*ct3 + c_m2(2,1)*beta_mat2;
y2 = c_m2(1,1)*ct3 + c_m2(2,1)*beta_mat2 + c_m2(3,1)*cT3 + c_m2(4,1)*cW3;
q2(1500,1)   % q1500
y2(1500,1)   % strain at t1500
%plot (y2,ct3)
ot = [el;el2] * (1/144) ;
ot_index = ot - ones(3,1); %outlier index

%% g bit

ac = [ct3 beta_mat2 cW3];
c_m3 = (ac'*ac)^-1*ac' * cs3 %coeff matrix

q3 = c_m3(1,1)*ct3 + c_m3(2,1)*beta_mat2;
y3 = c_m3(1,1)*ct3 + c_m3(2,1)*beta_mat2 + c_m3(3,1)*cW3;

q3(1500,1) %q1500
