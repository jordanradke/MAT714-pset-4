
pp = 6;
err = zeros(pp,1);
h   = zeros(pp,1);

Nbench = 2^10-1;
dxbench = 1/Nbench;
dtbench = .38*1/Nbench;
XXX = 0:dxbench:1;
TTT  = 0:dtbench:10;
ubench = godunov(Nbench,10);
[Tgrid,Xgrid] = meshgrid(TTT,XXX);

for p = 2:pp
    p
    % grid refinement + time steps
    N = 2^p-1;
    T = 10;
    dx = 1/N;
    dt = .38/N;

    x = 0:dx:1;
    t = 0:dt:10;
    [TT,XX] = meshgrid(t,x);

    % numerical solution
    fprintf('numerical solution...')
    u = burger(N,T);
    
    % benchmark solution sampled on grid (only for f=0 init. data)
    fprintf('benchmark solution...')
    u_sol = interp2(Tgrid,Xgrid,ubench,TT,XX);
    
    err(p) = max(max(max(abs(u-u_sol))));
    h(p) = dx;
    
end

% plot loglog error 
figure(1); clf();
loglog(h,err,'o-', 'LineWidth', 2)
hold on;
loglog(h, h.^2, 'LineStyle', '-')

ax = gca;
ax.YAxis.FontSize = 13;
ax.XAxis.FontSize = 13;

title('Error','FontSize',24);
xlabel('$h$','Interpreter','latex','FontSize',24)
ylabel('relative $\ell^\infty$ error', 'Interpreter','latex','FontSize',24)
    
    

    
    