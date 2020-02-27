clc,clear;
P=50;%total number of points (workers)
N=47;%interpolation degree+1 (points needed)

p=1:P;
grid1=2*(p-P/2)./P; %equi-distant grid
grid2=cos(pi*((2*p-1))./(2*P)); %chebyshev grid


if mod(P,2)
pp=0:(P-1)/2-1;
gg=1+pp./2;
grid3=[-gg(end:-1:1) 0 gg];
else
pp=0:P/2-1;   
gg=1+pp./2;
grid3=[-gg(end:-1:1)  gg];
end



M_g1=zeros(P,N);
M_g2=zeros(P,N);
A_g1=zeros(P,N);
A_g2=zeros(P,N);
M_g3=zeros(P,N);

%evaluation matrices w/ different grids
for i=1:N
M_g1(:,i)=(grid1.^(i-1))';
M_g2(:,i)=(grid2.^(i-1))';
M_g3(:,i)=(grid3.^(i-1))';
A_g1(:,i)=(cos((i-1).*acos(grid1)))';
A_g2(:,i)=(cos((i-1).*acos(grid2)))';
end

choices=nchoosek(p,N);

cond_M_g1_max=cond(M_g1(choices(1,:),:));
cond_M_g2_max=cond(M_g2(choices(1,:),:));
cond_A_g1_max=cond(A_g1(choices(1,:),:));
cond_A_g2_max=cond(A_g2(choices(1,:),:));
cond_M_g3_max=cond(M_g3(choices(1,:),:));

cond_M_g1_avg=cond_M_g1_max/length(choices);
cond_M_g2_avg=cond_M_g2_max/length(choices);
cond_A_g1_avg=cond_A_g1_max/length(choices);
cond_A_g2_avg=cond_A_g2_max/length(choices);
cond_M_g3_avg=cond_M_g3_max/length(choices);

for i=2:length(choices)
    i  
current_cond_M_g1=cond(M_g1(choices(i,:),:));   
current_cond_M_g2=cond(M_g2(choices(i,:),:));   
current_cond_A_g1=cond(A_g1(choices(i,:),:));   
current_cond_A_g2=cond(A_g2(choices(i,:),:));   
current_cond_M_g3=cond(M_g3(choices(i,:),:));   
 
cond_M_g1_avg= cond_M_g1_avg+current_cond_M_g1/length(choices);
cond_M_g2_avg= cond_M_g2_avg+current_cond_M_g2/length(choices);
cond_A_g1_avg= cond_A_g1_avg+current_cond_A_g1/length(choices);
cond_A_g2_avg= cond_A_g2_avg+current_cond_A_g2/length(choices);
cond_M_g3_avg= cond_M_g3_avg+current_cond_M_g3/length(choices);



if current_cond_M_g1> cond_M_g1_max
cond_M_g1_max=cond(M_g1(choices(i,:),:));
end
if current_cond_M_g2> cond_M_g2_max
cond_M_g2_max=cond(M_g2(choices(i,:),:));

end
if current_cond_A_g1> cond_A_g1_max
cond_A_g1_max=cond(A_g1(choices(i,:),:));
end
if current_cond_A_g2> cond_A_g2_max
cond_A_g2_max=cond(A_g2(choices(i,:),:));
end

if current_cond_M_g3> cond_M_g3_max
cond_M_g3_max=cond(M_g3(choices(i,:),:));
end
end

display(cond_A_g2_max),display(cond_A_g2_avg)
display(cond_A_g1_max),display(cond_A_g1_avg)
display(cond_M_g2_max),display(cond_M_g2_avg)
display(cond_M_g1_max),display(cond_M_g1_avg)
display(cond_M_g3_max),display(cond_M_g3_avg)
