%% This function finds the 2nd Order Finite Difference of the LD Field

function [DLD] = LG_FD_EXO2(LD,dx,dy)
[ny, nx] = size(LD);
Dx = zeros(ny,nx);
Dy = zeros(ny,nx);

for i = 1:ny
    Dx(i,:) = FD_EXO2(LD(i,:), dx);
end
for i = 1:nx
    Dy(:,i) = FD_EXO2(LD(:,i), dy);
end
DLD.LDX = Dx;
DLD.LDY = Dy;