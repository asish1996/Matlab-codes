nel=500;

nnp= 2*nel + 1;
K=zeros(nnp,nnp);
Kd=zeros(nnp,nnp)
fl=zeros(nnp,1);
x=zeros(nnp,1);
Lel=1/nel;

%stiffness matrices
for i=3:2:nnp
    range=[i-2 i-1 i];
    K(range,range)=K(range,range)+(1/(3*Lel))*[7 -8 1;-8 16 -8;1 -8 7];
end

for i=3:2:nnp
    range=[i-2 i-1 i];
    Kd(range,range)=Kd(range,range)+(Lel/2)*[0.2667 0.1333 -0.0667;0.1333 1.0667 0.1333;-0.0667 0.1333 0.2667];
end


Kf = K - Kd

%Boundary vector
a=[1];
b=zeros(nnp-1,1);
fb = [b;a];

%partitioning of matrices
Kt=Kf(2:nnp,2:nnp);
fbt=fb(2:nnp);
s = inv(Kt)*fbt;

s=[0;s];

% X-axis
for i = 1 : nnp
    range = [i]
    x(range,1)=x(range,1)+ [(i-1)*(Lel/2)];
end

x

plot(x,s)
hold on

fplot(@(x) 1.851*sin(x),[0 1],'r')







