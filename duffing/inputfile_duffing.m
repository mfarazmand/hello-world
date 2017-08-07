clc; 
clear all;

m=50; n=100;
x = linspace(-2,2,m);
y = linspace(-1,1,n);
[xi, yi]=meshgrid(x,y);
x0=reshape(xi,[],1);
y0=reshape(yi,[],1);
v1_1=ones(length(x0),1);
v1_2=zeros(length(x0),1);

% parameters
T=20;
tspan = linspace(0,T,3);
X0=[x0;y0;v1_1;v1_2];

% rho=1e-5;
% [F11, F12, F21, F22, l1, l2, v1, v2]=cg_strain_tensor(@duffing , xi, yi, tspan, rho);
% [C11, C12, C22, l1, l2, v1, v2]=DF2C(F11,F12,F21,F22);


options = odeset('RelTol',1e-5,'AbsTol',1e-5);
[time,F]=ode45(@duffing_do,tspan,X0,options);

N=length(x0);
xt=F(end,1:N)'; yt=F(end,N+1:2*N)';
v1_1=F(end,2*N+1:3*N)'; v1_2=F(end,3*N+1:4*N)';
v1_1=reshape(v1_1,n,m);
v1_2=reshape(v1_2,n,m);
quiver(xi,yi,v1_1,v1_2,'k')