clc,clear;
P=67;%total number of points (workers)
R=P-3;%interpolation degree+1 (points needed)
m=sqrt(R);
s=2;% A_i,B_i have dimensions s x N2, and N2 x s.
N1=m*s;
N2=1000;
N3=N1;
p=1:P;

H=zeros(m^2,m^2);
H(1:m,1:m)=eye(m);
for i = 1:m-1
    H(i*m+1,i*m+1)=1;
    for j=1:m-1
    H(i*m+1+j,i*m+1+j)=1/2;
    H(i*m+1-j,i*m+1+j)=1/2;
    end
end

grid2=cos(pi*((2*p-1))./(2*P)); %chebyshev grid
grid2=grid2(end:-1:1);

A=normrnd(0,1,N1,N2);
B=normrnd(0,1,N2,N3);

M_g2=zeros(P,R);
C_g2=zeros(P,R);

Out_C=zeros(N1,N3);
Out_M=zeros(N1,N3);

%Evaluation Matrices for Mon. and Cheb.
for i=1:R
M_g2(:,i)=(grid2.^(i-1))';
C_g2(:,i)=(cos((i-1).*acos(grid2)))';
end
%C_g2(:,1)=0.5*C_g2(:,1);

pA_C=zeros(s,N2,P); 
pB_C=zeros(N2,s,P);
pAB_C=zeros(s,s,P);

pA_M=zeros(s,N2,P); 
pB_M=zeros(N2,s,P);
pAB_M=zeros(s,s,P);



for i=1:P
for j=1:m
pA_C(:,:,i)=pA_C(:,:,i)+A((j-1)*s+1:j*s,:).*C_g2(i,j);
pB_C(:,:,i)=pB_C(:,:,i)+B(:,(j-1)*s+1:j*s).*C_g2(i,(j-1)*m+1);

pA_M(:,:,i)=pA_M(:,:,i)+A((j-1)*s+1:j*s,:).*M_g2(i,j);
pB_M(:,:,i)=pB_M(:,:,i)+B(:,(j-1)*s+1:j*s).*M_g2(i,(j-1)*m+1);
end

pAB_C(:,:,i)=pA_C(:,:,i)*pB_C(:,:,i); 
pAB_M(:,:,i)=pA_M(:,:,i)*pB_M(:,:,i); 
end

Coef_C_g2=zeros(s,s,R);
out=zeros(s,s,R);
Coef_M_g2=zeros(s,s,R);

choices=nchoosek(p,R);
inv_C_g2_max=inv(C_g2(choices(1,:),:));
inv_M_g2_max=inv(M_g2(choices(1,:),:));


for i=1:s
for j=1:s
Coef_C_g2(i,j,:)=inv_C_g2_max*squeeze(pAB_C(i,j,choices(1,:)));
out(i,j,:)=inv(H)*squeeze(Coef_C_g2(i,j,:));
Coef_M_g2(i,j,:)=inv_M_g2_max*squeeze(pAB_M(i,j,choices(1,:)));
end
end

%Concatenation of coefs.
for i=1:s
for j=1:s
Sh_C=reshape(out(i,j,:),m,m);
Sh_M=reshape(Coef_M_g2(i,j,:),m,m);
Out_C(i:s:N1-s+i,j:s:N1-s+j)=Sh_C;
Out_M(i:s:N1-s+i,j:s:N1-s+j)=Sh_M;
end
end

error_C_g2_max=sqrt(sum(sum((Out_C-A*B).^2)))/sqrt(sum(sum((A*B).^2)));
error_M_g2_max=sqrt(sum(sum((Out_M-A*B).^2)))/sqrt(sum(sum((A*B).^2)));


error_C_g2_avg=error_C_g2_max/nchoosek(P,R);
error_M_g2_avg=error_M_g2_max/nchoosek(P,R);


for k=2:length(choices)
    k
inv_C_g2_current=inv(C_g2(choices(k,:),:));
inv_M_g2_current=inv(M_g2(choices(k,:),:));


for i=1:s
for j=1:s
Coef_C_g2(i,j,:)=inv_C_g2_current*squeeze(pAB_C(i,j,choices(k,:)));
out(i,j,:)=inv(H)*squeeze(Coef_C_g2(i,j,:));
Coef_M_g2(i,j,:)=inv_M_g2_current*squeeze(pAB_M(i,j,choices(k,:)));
end
end

%Concatenation of coefs.
for i=1:s
for j=1:s
Sh_C=reshape(out(i,j,:),m,m);
Sh_M=reshape(Coef_M_g2(i,j,:),m,m);
Out_C(i:s:N1-s+i,j:s:N1-s+j)=Sh_C;
Out_M(i:s:N1-s+i,j:s:N1-s+j)=Sh_M;
end
end

error_C_g2_current=sqrt(sum(sum((Out_C-A*B).^2)))/sqrt(sum(sum((A*B).^2)));
error_M_g2_current=sqrt(sum(sum((Out_M-A*B).^2)))/sqrt(sum(sum((A*B).^2)));

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


