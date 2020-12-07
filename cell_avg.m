function [cells] = cell_avg(N)

M = 2^12;

dy = 1/(M-1);
y = 0:dy:1;


u = 1.5+sin(2*pi*y);
u = circshift(u, floor(M/(2*(N+1))));

% fit to array of size N+1, each averaged around n.
umat = reshape(u, [M/(N+1), N+1]);

cells = mean(umat);

end
