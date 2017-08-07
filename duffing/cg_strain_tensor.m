function [F11, F12, F21, F22, l1, l2, v1, v2]=cg_strain_tensor(rhs , xi, yi, tspan, rho);
%% COMMENTS
% rhs    = a function which represents the right hand side of the ode to be integrated
% xi, yi = the x and y coordinate of the initial position of the Lagrangian
%          particles
% tspan  = time span over which to integrate
% rho    = size of the local grid
%%
[m,n]=size(xi);
Nrad = 5;
xt=zeros(m,n,Nrad); yt=zeros(m,n,Nrad);
for k=1:Nrad-1
    xt(:,:,k) = xi + rho*cos((k-1)*pi/2);
    yt(:,:,k) = yi + rho*sin((k-1)*pi/2);
end
xt(:,:,Nrad) = xi; yt(:,:,Nrad)=yi;
%%
%%*************Time integration***********************************
xt_vec = reshape(xt, Nrad*m*n, 1);
yt_vec = reshape(yt, Nrad*m*n, 1);
options = odeset('RelTol',1e-7,'AbsTol',1e-7);
[~, F] = ode45(rhs, tspan, [xt_vec; yt_vec],options);
%%
[sizeT, N] = size(F);
xt = reshape(F(sizeT,1:N/2)      , m, n, Nrad);
yt = reshape(F(sizeT,N/2+1:N)  , m, n, Nrad);

%%
%*********computation of eigen-values and eigen-vectors
% l1=zeros(m,n,1); l2 = zeros(m,n, 1);
v1 = zeros(m,n,2); v2 = zeros(m,n,2);
size(xt)
size(yt)
F11 = (xt(:,:,1)-xt(:,:,3))/(2*rho);
F12 = (xt(:,:,2)-xt(:,:,4))/(2*rho);
F21 = (yt(:,:,1)-yt(:,:,3))/(2*rho);
F22 = (yt(:,:,2)-yt(:,:,4))/(2*rho);

C11=F11.^2+F21.^2;
C12=F11.*F12+F22.*F21;
C22=F22.^2+F12.^2;

trC=C11+C22;
detC=C11.*C22-C12.^2;

l1 = 0.5*trC-sqrt((0.5*trC).^2-detC);
l2 = 0.5*trC+sqrt((0.5*trC).^2-detC);

v2(:,:,1) = -C12./sqrt(C12.^2+(C11-l2).^2);
v2(:,:,2) = (C11-l2)./sqrt(C12.^2+(C11-l2).^2);
v1(:,:,1) = v2(:,:,2);
v1(:,:,2) =-v2(:,:,1);

%% l1 and l2 on the regular grid
% xt_reg=zeros(m,n); yt_reg=zeros(m,n);
% xt_reg(:,:)=xt(:,:,Nrad); yt_reg(:,:)=yt(:,:,Nrad);
% [gradF11,gradF12]= gradient(xt_reg,xi',yi);
% [gradF21,gradF22]= gradient(yt_reg,xi',yi);
% 
% C11=gradF11.^2+gradF21.^2;
% C12=gradF11.*gradF12+gradF22.*gradF21;
% C22=gradF22.^2+gradF12.^2;
% 
% trC=C11+C22;
% detC=C11.*C22-C12.^2;
% 
% l1 = 0.5*trC-sqrt((0.5*trC).^2-detC);
% l2 = 0.5*trC+sqrt((0.5*trC).^2-detC);