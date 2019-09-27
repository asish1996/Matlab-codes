nel=500;

nel1=nel/2;
nnp1=2*nel1+1;
nnp=2*nel+1;
K=zeros(nnp1,nnp1);
k=zeros(nnp1,nnp1);
fl=zeros(nnp,1);
fb=zeros(nnp,1);
x=zeros(nnp,1);
Lel=(1/nel)*0.7;
Lel1=(1/nel1)*0.3;
Lel2=(1/nel1)*0.4;

%Stiffness matrix
for i=3:2:nnp1
    range=[i-2 i-1 i];
    K(range,range)=K(range,range)+(1/(3*Lel1))*[7 -8 1;-8 16 -8;1 -8 7];
end
K
for i=3:2:nnp1
    range=[i-2 i-1 i];
    k(range,range)=k(range,range)+(4/(3*Lel2))*[7 -8 1;-8 16 -8;1 -8 7];
end

k

%Load vector
for i=3:2:nnp
    range=[i-2;i-1;i];
    fl(range,1)=fl(range,1)+ ((Lel*Lel)/6)*[i-3;i-3;i-2];
end

fl

li=nnp;
c=zeros(nnp);
for i = size(k,1):-1:1
    lj = nnp;
    for j = size(k,1):-1:1
        c(li,lj)= k(i,j);
        lj = lj-1;
    end
    li = li- 1;

end
c

P=K;
for i = size(P,1)+1:nnp
    P(i,:)=0;
    P(:,i)=0;
end

P

KS=P+c;

KSt = KS(2:nnp-1,2:nnp-1);
flt = fl(2:nnp-1,1);

s = inv(KSt)*flt;
s(nnp-1,1) = 0;
r=[0];
ss=[r;s]

%X-axis
x1=zeros(2*nel1+1,1);
x2=zeros(2*nel1,1);

for i = 1 : 2*nel1+1
    range = [i]
    x1(range,1)=x1(range,1)+ [(i-1)*(Lel1/2)];
end

for i= 1 : 2*nel1 
    range = [i]
    x2(range,1)=x2(range,1)+ [x1(2*nel1+1,1)+(i)*(Lel2/2)];
end

x = [x1;x2];
x
% du/dx
% for i=1:2:nnp-1
%     p = (3/Lel)*(ss(i+1,1) - ss(i,1));
%     q = [p p];
%     r = [x(i,1) x(i+2,1)];
%     plot (r,q)
%     hold on
% end

 plot(x,ss)
 hold on





    

