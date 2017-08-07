function [C11, C12, C22, l1, l2, v1, v2]=DF2C(F11,F12,F21,F22)

C11=F11.^2+F21.^2;
C12=F11.*F12+F22.*F21;
C22=F22.^2+F12.^2;

trC=C11+C22;
detC=C11.*C22-C12.^2;

l2 = 0.5*trC+sqrt((0.5*trC).^2-detC);
% l1 = detC./( 0.5*trC+sqrt((0.5*trC).^2-detC) );
l1 = 0.5*trC-sqrt((0.5*trC).^2-detC);

v2(:,:,1) = -C12./sqrt(C12.^2+(C11-l2).^2);
v2(:,:,2) = (C11-l2)./sqrt(C12.^2+(C11-l2).^2);
v1(:,:,1) = v2(:,:,2);
v1(:,:,2) =-v2(:,:,1);