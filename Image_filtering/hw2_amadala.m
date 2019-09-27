Phi = double(imread('foot.pgm'));
%Phi = imread('foot.pgm');
A = imread('foot.pgm');
imwrite(A,'actual_image.pgm')
%imshow(A);

%% FTCS Explicit scheme

Phi_new = Phi;

dt = 0.1;
t = 2;

for k = 1:500;  %No.of time steps
    for i = 2:135;
        for j = 2:148;
            Phi_new(i,j) = Phi(i,j) + dt*(Phi(i+1,j)-2*Phi(i,j)+Phi(i-1,j)+Phi(i,j+1)-2*Phi(i,j)+Phi(i,j-1));
        end
    end
    Phi = Phi_new;
end

round_array = uint8(Phi_new);
%imshow(round_array)
imwrite(round_array,'linearfilter_image.pgm')


%% Non linear filtering

Phi2 = double(A);
Phi2_new = Phi2;
lambda = 10;
dt = 0.01;
%t = 1;

for k = 1:500;      %No.of time steps
    for i = 2:135;
        for j = 2:148;
            g_delphi = abs((Phi2(i+1,j)-Phi2(i-1,j)+Phi2(i,j+1)-Phi2(i,j-1))/2);
            Phi2_new(i,j) = Phi2(i,j) + dt*(Phi2(i-1,j)+Phi2(i+1,j)+Phi2(i,j-1)+Phi2(i,j+1)-4*Phi2(i,j))/(1+((g_delphi^2)/lambda^2));       
        end
    end
    Phi2 = Phi2_new;
end

round_array2 = uint8(Phi2_new);
imshow(round_array2)
imwrite(round_array2,'nonlinearfilter_image.pgm')
            
