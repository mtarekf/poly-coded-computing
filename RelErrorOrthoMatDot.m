clc,clear;
P=50;%total number of points (workers)
R=49;%interpolation degree+1 (points needed)
N1=10;
m=(R+1)/2;
s=2;% A_i,B_i have dimensions N1 x s, and s x N3.
N2=m*s;
N3=10;
p=1:P;
grid2=cos(pi*((2*p-1))./(2*P)); %chebyshev grid
grid2=grid2(end:-1:1);

A=normrnd(0,1,N1,N2);
B=normrnd(0,1,N2,N3);

M_g2=zeros(P,R);
C_g2=zeros(P,R);


%Evaluation Matrices for Mononmial basis and Chebyshev basis
for i=1:R
M_g2(:,i)=(grid2.^(i-1))';
C_g2(:,i)=(cos((i-1).*acos(grid2)))';
end
C_g2(:,1)=0.5*C_g2(:,1);

pA_C=zeros(N1,s,P); 
pB_C=zeros(s,N3,P);
pAB_C=zeros(N1,N3,P);

pA_M=zeros(N1,s,P); 
pB_M=zeros(s,N3,P);
pAB_M=zeros(N1,N3,P);



for i=1:P

for j=1:m
pA_C(:,:,i)=pA_C(:,:,i)+A(:,(j-1)*s+1:j*s).*C_g2(i,j);
pB_C(:,:,i)=pB_C(:,:,i)+B((j-1)*s+1:j*s,:).*C_g2(i,m-j+1);

pA_M(:,:,i)=pA_M(:,:,i)+A(:,(j-1)*s+1:j*s).*M_g2(i,j);
pB_M(:,:,i)=pB_M(:,:,i)+B((j-1)*s+1:j*s,:).*M_g2(i,m-j+1);
end

pAB_C(:,:,i)=pA_C(:,:,i)*pB_C(:,:,i); 
pAB_M(:,:,i)=pA_M(:,:,i)*pB_M(:,:,i); 
end

Coef_C_g2=zeros(N1,N3,R);
Coef_M_g2=zeros(N1,N3,R);

choices=nchoosek(p,R);
inv_C_g2_max=inv(C_g2(choices(1,:),:));
inv_M_g2_max=inv(M_g2(choices(1,:),:));


for i=1:N1
for j=1:N3
Coef_C_g2(i,j,:)=2*(inv_C_g2_max*squeeze(pAB_C(i,j,choices(1,:))));
Coef_M_g2(i,j,:)=inv_M_g2_max*squeeze(pAB_M(i,j,choices(1,:)));
end
end
error_C_g2_max=sqrt(sum(sum((Coef_C_g2(:,:,m)-A*B).^2)))/sqrt(sum(sum((A*B).^2)));
error_M_g2_max=sqrt(sum(sum((Coef_M_g2(:,:,m)-A*B).^2)))/sqrt(sum(sum((A*B).^2)));


error_C_g2_avg=error_C_g2_max/nchoosek(P,R);
error_M_g2_avg=error_M_g2_max/nchoosek(P,R);

for k=2:length(choices)
    k
inv_C_g2_current=inv(C_g2(choices(k,:),:));
inv_M_g2_current=inv(M_g2(choices(k,:),:));


for i=1:N1
for j=1:N3
Coef_C_g2(i,j,:)=2*(inv_C_g2_current*squeeze(pAB_C(i,j,choices(k,:))));
Coef_M_g2(i,j,:)=inv_M_g2_current*squeeze(pAB_M(i,j,choices(k,:)));
end
end
error_C_g2_current=sqrt(sum(sum((Coef_C_g2(:,:,m)-A*B).^2)))/sqrt(sum(sum((A*B).^2)));
error_M_g2_current=sqrt(sum(sum((Coef_M_g2(:,:,m)-A*B).^2)))/sqrt(sum(sum((A*B).^2)));

error_C_g2_avg=error_C_g2_avg+error_C_g2_current/nchoosek(P,R);
error_M_g2_avg=error_M_g2_avg+error_M_g2_current/nchoosek(P,R);

if error_C_g2_current> error_C_g2_max
error_C_g2_max=error_C_g2_current;
end
if error_M_g2_current> error_M_g2_max
error_M_g2_max=error_M_g2_current;
end
end

error_C_g2_max
error_M_g2_max
error_C_g2_avg
error_M_g2_avg

