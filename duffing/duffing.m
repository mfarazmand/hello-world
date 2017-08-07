function dy = duffing(t,y)
t
delta=0.00;
tau=0.00;
N=round(length(y)/2);
dy = zeros(2*N,1);    % a column vector
dy(1:N,1) = y(N+1:2*N,1);
dy(N+1:2*N,1) = (-delta*y(N+1:2*N,1)+y(1:N,1)-y(1:N,1).^3+tau*cos(t));


% dy(1:N,1) = y(N+1:2*N,1);
% dy(N+1:2*N,1) = -sin(y(1:N,1));