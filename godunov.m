%% godunov method for burger's equation
function [w] = godunov(N,T)

%N = 19;
%T = 10;

dx = 1/N;
dt = .38/N;  % CFL stability condition: dt < 2/5*dx 
TT = ceil(T/dt);

q = zeros(N+1,TT);

% initial conditions @ t=0
q(:,1) = cell_avg(N);
%q(:,1) = 1.5 + sin(2*pi*x);

for t = 2:TT
    % enforce periodic boundary conditions
    q(1,t)   = q(1,t-1) - dt/dx*(flux(q(1,t-1),q(2,t-1)) - flux(q(N+1,t-1),q(1,t-1)));
    q(N+1,t) = q(N+1,t-1) - dt/dx*(flux(q(N+1,t-1),q(1,t-1)) - flux(q(N,t-1),q(N+1,t-1)));

    for i = 2:N
        % flux differencing formula to update cell averages
        q(i,t) = q(i,t-1) - dt/dx*(flux(q(i,t-1),q(i+1,t-1))-flux(q(i-1,t-1),q(i,t-1)));
    end  
end

%surf(1:floor(TT),x,q(:,1:floor(TT)))

w = q;

end

