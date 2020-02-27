clc,clear;
P=20;%total number of points (workers)
N=18;%interpolation degree+1 (points needed)
d=10;
%m=(N+1)/2;
p=1:P;
grid2=cos(pi*((2*p-1))./(2*P)); %chebyshev grid
grid2=grid2(end:-1:1);
A=normrnd(0,1,d,N);
B=normrnd(0,1,d,1);

M_g2=zeros(P,N);
C_g2=zeros(P,N);
%Evaluation Matrices for Mon. and Cheb.
for i=1:N
M_g2(:,i)=(grid2.^(i-1))';
C_g2(:,i)=(cos((i-1).*acos(grid2)))';
end

%Polynomial Evaluations for Lagrange
L=ones(N,P-N);
for k=N+1:P
for i=1:N    
Denominator=1;
Numerator=1;
for j=1:N
if j==i
    %Denominator=Denominator;
    %Numerator=Numerator;
else
    Denominator=Denominator*(grid2(1,i)-grid2(1,j));
    Numerator=Numerator*(grid2(1,k)-grid2(1,j));
end
end  
L(i,k)=Numerator/Denominator;
end
end
L_g2A=zeros(d,P);
L_g2A(:,1:N)=A(:,1:N);
for k=N+1:P
    for i=1:N
L_g2A(:,k)=L_g2A(:,k)+A(:,i)*L(i,k);
    end
end

f_L_g2=B'*L_g2A;


choices=nchoosek(p,N);
inv_M_g2_max=inv(M_g2(choices(1,:),:));
inv_C_g2_max=inv(C_g2(choices(1,:),:));

Coef_M_g2=inv_M_g2_max*f_L_g2(:,choices(1,:))';
Coef_C_g2=inv_C_g2_max*f_L_g2(:,choices(1,:))';

Out_M_g2=M_g2(1:N,:)*Coef_M_g2;
Out_C_g2=C_g2(1:N,:)*Coef_C_g2;


error_M_g2_max=sqrt(sum((Out_M_g2'-f_L_g2(1:N)).^2))/sqrt(sum((f_L_g2(1:N)).^2));
error_C_g2_max=sqrt(sum((Out_C_g2'-f_L_g2(1:N)).^2))/sqrt(sum((f_L_g2(1:N)).^2));


error_M_g2_avg=error_M_g2_max/nchoosek(P,N);
error_C_g2_avg=error_C_g2_max/nchoosek(P,N);

for i=2:length(choices)
    i
inv_M_g2_current=inv(M_g2(choices(i,:),:));
inv_C_g2_current=inv(C_g2(choices(i,:),:));

Coef_M_g2=inv_M_g2_current*f_L_g2(:,choices(i,:))';
Coef_C_g2=inv_C_g2_current*f_L_g2(:,choices(i,:))';

Out_M_g2=M_g2(1:N,:)*Coef_M_g2;
Out_C_g2=C_g2(1:N,:)*Coef_C_g2;


error_M_g2_current=sqrt(sum((Out_M_g2'-f_L_g2(1:N)).^2))/sqrt(sum((f_L_g2(1:N)).^2));
error_C_g2_current=sqrt(sum((Out_C_g2'-f_L_g2(1:N)).^2))/sqrt(sum((f_L_g2(1:N)).^2));

error_M_g2_avg=error_M_g2_avg+error_M_g2_current/nchoosek(P,N);
error_C_g2_avg=error_C_g2_avg+error_C_g2_current/nchoosek(P,N);

if error_M_g2_current> error_M_g2_max
error_M_g2_max=error_M_g2_current;
end

if error_C_g2_current> error_C_g2_max
error_C_g2_max=error_C_g2_current;
end
end

display(error_C_g2_max)
display(error_M_g2_max)
display(error_C_g2_avg)
display(error_M_g2_avg)


