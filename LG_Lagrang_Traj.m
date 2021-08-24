% To advect a defined grid inside the instantenous velocity field and to
% generate the Lagrangian trajectory​
function [A_FWD,A_BWD,I,VECRF,nx,ny] = LG_Lagrang_Traj(VEC,t,SRF,indx,tau_fwd,tau_bwd,varargin)


Szx = size(VEC{1}.X,2);
x = linspace(VEC{1}.X(1),VEC{1}.X(end),SRF*Szx);

Szy = size(VEC{1}.Y,1);
y = linspace(VEC{1}.Y(1),VEC{1}.Y(end),SRF*Szy);

[P.X,P.Y] = meshgrid(x,y);
VECRF.C = logical(interp2(VEC{1}.X, VEC{1}.Y, double(VEC{1}.C), P.X, P.Y));
VECRF.X = interp2(VEC{1}.X, VEC{1}.Y, VEC{1}.X, P.X, P.Y);
VECRF.Y = interp2(VEC{1}.X, VEC{1}.Y, VEC{1}.Y, P.X, P.Y);
% I = find(VEC{1}.C ==1);
% k = sum(VEC{1}.C,'all');
I = find(VECRF.C ==1);
k = sum(VECRF.C,'all');
G.X = VECRF.X(I);
G.Y = VECRF.Y(I);
% count =1;
% while count <= k
%     %     G.X(count) = VEC{1}.X(I(count));
%     %     G.Y(count) = VEC{1}.Y(I(count));
%     G.X(count) = VECRF.X(I(count));
%     G.Y(count) = VECRF.Y(I(count));
%     count = count+1;
% end
% nx = size(VEC{1}.X, 2);
% ny = size(VEC{1}.Y, 1);
nx = size(VECRF.X, 2);
ny = size(VECRF.Y, 1);
% addpath('C:\Users\a_arwis\Ahmed PhD\MATfluids-master-Nov') % To include the
% directory to PP_Advectpoints from MATfluids package
if nargin == 7 && strcmpi(varargin{1},'append')
    VEC = [VEC; VEC];
    G_FWD = PP_AdvectGrid(G, VEC(indx:indx+tau_fwd), t(2)-t(1),'singlestep');
    G_BWD = PP_AdvectGrid(G, flipud(VEC(indx-tau_bwd+(0.5*length(VEC)):indx+(0.5*length(VEC)))), t(1)-t(2),'singlestep');
else
    G_FWD = PP_AdvectGrid(G, VEC(indx:indx+tau_fwd), t(2)-t(1),'singlestep');
    G_BWD = PP_AdvectGrid(G, flipud(VEC(indx-tau_bwd:indx)), t(1)-t(2),'singlestep');
end
% Organize the output in a 2D matrix form for LD computations
A_BWD.X = extractField(G_BWD,'X')';
A_BWD.Y = extractField(G_BWD,'Y')';
A_FWD.X = extractField(G_FWD,'X')';
A_FWD.Y = extractField(G_FWD,'Y')';
