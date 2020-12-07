%% upwind method for inviscid burger's eqn
function [w] = burger(N,T)

dx = 1/N;
dt = .38/N;  % CFL stability condition: dt < 2/5*dx 
TT = ceil(T/dt);


x = 0:dx:1;
u = zeros(N+1,TT);

% initial conditions @ t=0
u(:,1) = 1.5 + sin(2*pi*x);


for t = 2:floor(T/dt)    
    % enforce periodic boundary conditions
    fluxL = max(u(1,t-1),0).*(u(1,t-1) - u(N+1,t-1)) ...
            + min(u(1,t-1),0).*(u(2,t-1) - u(1,t-1));
    u(1,t) = u(1,t-1) - dt/dx.*fluxL;
        
    fluxR = max(u(N+1,t-1),0).*(u(N+1,t-1) - u(N,t-1)) ...
            + min(u(N+1,t-1),0).*(u(1,t-1) - u(N+1,t-1));
    u(N+1,t) = u(N+1,t-1) - dt/dx.*fluxR;
    
    % main interior loop
    flux = max(u(2:N,t-1),0).*(u(2:N,t-1) - u(1:N-1,t-1)) ... 
           + min(u(2:N,t-1),0).*(u(3:N+1,t-1) - u(2:N,t-1));
    u(2:N,t) = u(2:N,t-1) - dt/dx.*flux;
end

w = u;

end
