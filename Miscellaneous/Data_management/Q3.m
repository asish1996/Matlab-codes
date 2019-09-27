load Water_Temp_Strain.txt
%load cleaned_data1.csv

sensors = Water_Temp_Strain(:,:);
W = sensors(:,1); %Water level in meters
T = sensors(:,2); %Temperature in celsius
S = sensors(:,3); %Strain in micro starins
S = S %* 10^-6;

nrows = [1:3650];
rowindex = [0;nrows'];
delT = 144;
time = rowindex.*delT;

%% p=1 & n=3;

sn1 = zeros(3649,1);
sn = [1.677454715388890/100000;sn1;4.807054578/1000000];

for i=2:3650;
    
    sn(i,1) = (S(i-1,1) + S(i,1) + S(i+1,1))/3;
end

px1 = smoothts(sn','b',3)
%plot(time,px1')
%plot (time,sn)

%% p=4 & n=9

sn2 = [S(1:5,1);zeros(3641,1);S(3647:3651,1)];

for i=5:3646;
    sn2(i,1) = (S(i-4,1)+S(i-3,1)+S(i-2,1)+S(i-1,1)+S(i,1)+S(i+1,1)...
        +S(i+2,1)+S(i+3,1)+S(i+4,1))/9;
end
%plot(time,sn2)
px2 = smoothts(sn2','b',9);
%plot(time,px2')
%% p=20 & n=41

sn3 = S;

for i=21:3630;
    sn3(i,1) = sum(sn3(i-20:i+20,1))/41;
end

%plot(time,sn3)
px3 = smoothts(sn3','b',41);
%plot(time,px3')

% sn(1500,1)
% sn2(1500,1)
% sn3(1500,1)

% nrows = [1:3650];
% rowindex = [0;nrows'];
% delT = 144;
% time = rowindex.*delT; 
% beta_mat = ones(3651);
% beta_mat = beta_mat(:,1);
% 
% A = [time beta_mat T W];
% c_m1 = (A'*A)^-1*A' * sn;
% q1 = c_m1(1,1)*time + c_m1(2,1)*beta_mat;
% %y = c_m1(1,1)*time + c_m1(2,1)*beta_mat + c_m1(3,1)*T + c_m1(4,1)*W;
% c_m2 = (A'*A)^-1*A' * sn2;
% q2 = c_m2(1,1)*time + c_m2(2,1)*beta_mat;
% c_m3 = (A'*A)^-1*A' * sn3;
% q3 = c_m3(1,1)*time + c_m3(2,1)*beta_mat;
% 
% q1(1500,1)
% q2(1500,1)
% q3(1500,1)

%% Smoothening using cleaned data

 X = table2array(cleaneddata1);
 tc = X(:,2);
 st = X(:,1);
 Tc = X(:,3);
 Wc = X(:,4);
 
 %% p=1 & n=3;

snc1 = st;

for i=2:3647;    
    snc1(i,1) = (st(i-1,1) + st(i,1) + st(i+1,1))/3;
end

snc1;
px4 = smoothts(st','b',3)
%plot(tc,st)
%plot (tc,px4')
%% p=4 & n=9

snc2 = st;

for i=5:3643;
    snc2(i,1) = (st(i-4,1)+st(i-3,1)+st(i-2,1)+st(i-1,1)+st(i,1)+st(i+1,1)...
        +st(i+2,1)+st(i+3,1)+st(i+4,1))/9;
end
%%plot(tc,snc2)
px5 = smoothts(st','b',9);
%plot(tc,px5')

%% p=20 & n=41

snc3 = st;

for i=21:3627;
    snc3(i,1) = sum(st(i-20:i+20,1))/41;
end

%plot(tc,snc3)
px6 = smoothts(st','b',41);
%plot(tc,px6')

% Aa = [tc beta_mat(1:3648,1) Tc Wc];
% c_m11 = (Aa'*Aa)^-1*Aa' * px4';
% q11 = c_m11(1,1)*tc + c_m11(2,1)*beta_mat(1:3648,1);
% %y = c_m1(1,1)*time + c_m1(2,1)*beta_mat + c_m1(3,1)*T + c_m1(4,1)*W;
% c_m22 = (Aa'*Aa)^-1*Aa' * px5';
% q22 = c_m22(1,1)*tc + c_m22(2,1)*beta_mat(1:3648,1);
% c_m33 = (Aa'*Aa)^-1*Aa' * px6';
% q33 = c_m33(1,1)*tc + c_m33(2,1)*beta_mat(1:3648,1);



