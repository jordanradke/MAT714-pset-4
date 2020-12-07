%% analytic solution of heat eq. via fourier on [0,1].
% This is not producing smooth solutions.

M = 20;
N = 10000;
T = 10;
dx = 1/N;
dt = 1/N;

nmodes = M;
a = zeros(nmodes,1);
b = zeros(nmodes,1);

% the function
e = .1;
xfun = @(x) exp(-1/(2*e)*(3.*x/2 + 1/(2*pi)*(1-cos(2*pi.*x))));

% numerically finding fourier coefficients
a0 = integral(xfun,0,1);

for m = 1:nmodes
    afun = @(x) exp(-1/(2*e)*(3.*x/2 + 1/(2*pi)*(1-cos(2*pi.*x)))).*cos(2*m*x.*pi);
    bfun = @(x) exp(-1/(2*e)*(3.*x/2 + 1/(2*pi)*(1-cos(2*pi.*x)))).*sin(2*m*x.*pi);
        
    a(m) = 2*integral(afun,0,1);
    b(m) = 2*integral(bfun,0,1);
end

% series solution for u
u = zeros(N+1,T);

u(:,1) = xfun(0:dx:1);

% is there a way to write in this via matrix operations?
for t = 2:T+1
    for i = 1:N+1
            % sum fourier series
            sum = a0;
            for m = 1:nmodes
                    sum = sum + a(m)*exp(-4*pi^2*m^2*(dt*(t-1)))*cos(2*pi*m*(dx*(i-1))) ...
                              + b(m)*exp(-4*pi^2*m^2*(dt*(t-1)))*sin(2*pi*m*(dx*(i-1)));
            end
            u(i,t) = sum;
    end
end

surf(0:1/T:1,0:1/N:1,u)

w = u;

% compute numerical derivative of heat eq solution
du = zeros(size(u));

for t = 1:T+1
    du(1,t)  = 1/(2*N)*(-3*u(1,t)+2*u(2,t) - u(3,t));
    du(N+1,t) = 1/(2*N)*(3*u(N-1,t) - 2*u(N,t) + u(N+1,t));
    
    du(2:N,t) = 1/(2*N)*(u(3:N+1,t)-u(1:N-1,t));
end

% hopf-cole transformation back 
surf(0:1/T:1,0:1/N:1, -2*e*du./u)

                    
 
