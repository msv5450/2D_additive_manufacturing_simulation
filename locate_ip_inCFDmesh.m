function NALU_eID = locate_ip_inCFDmesh(xy_ip, h, nx, ny)

% locates the NALU element that the given integration point belongs to
x_ip = xy_ip(1);
y_ip = xy_ip(2);

x_ip = x_ip + 0.5*nx*h;
y_ip = y_ip + 0.5*ny*h;
n_i = floor(x_ip/h);
n_j = floor(y_ip/h);
NALU_eID = n_j*nx + n_i + 1;

end