function dy = duffing_do(t,y)
t
delta=0.00;
tau=0.00;
N=round(length(y)/4);
dy = zeros(4*N,1);    % a column vector

x1 = y(1:N,1);
x2 = y(N+1:2*N,1);
v1_1 = y(2*N+1:3*N,1);
v1_2 = y(3*N+1:4*N,1);

dy(1:N,1) = x2;
dy(N+1:2*N,1) = x1-x1.^3+tau*cos(t);
dy(2*N+1:3*N,1) = v1_2-((2-3*x1.^2).*v1_1.*v1_2).*v1_1;
dy(3*N+1:4*N,1) = (1-3*x1.^2).*v1_1 - ((2-3*x1.^2).*v1_1.*v1_2).*v1_2;