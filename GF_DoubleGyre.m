%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                     Generate Flow : Double Gyre                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Laboratory of Cardiovascular Fluid Dyanmics                            %
%  Modified By :  Ahmed Darwish                                           %
%  Updated     :  16 April 2021                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          LIST OF VARIABLES                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INPUT ARGUMENTS 
%   tincr   : holds the time increment between t0 = 0 and T = 60          %
%   Ry      : holds the spatial resolution of the velocity grid (number of
%             points).
%%%% OUTPUT ARGUMENTS
%   x       : 2D grid of x coordinates
%   y       : 2D grid of y coordinates
%   t       : Vector of the simulated discrete time setps
%   VEC     : Cell of structures holding the instantenous velocity fields
%             and the spatial grid info and the mask (in VEC{#}.C).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x,y,t,VEC] = GF_DoubleGyre(tincr,Ry,varargin)
%%%%%%%%%%%% Variables and Constants %%%%%%%%%%%%%%%%
x    = linspace(0, 2, 2*Ry+1)';
y    = linspace(0, 1, Ry+1)';
t    = linspace(0,30,tincr)';
A    = 0.1;
epsn = 0.25;
omega = 0.5;

%%%%%%%%%%%% Generate the Flow Field %%%%%%%%%%%%%%
for k = 1:length(t)
    [VEC{k,1}.X,VEC{k,1}.Y] = meshgrid(x,flipud(y));
    a = epsn*sin(omega*t(k));
    b = 1 - 2*epsn*sin(omega*t(k));
    f = a*VEC{k}.X.^2 + b*VEC{k}.X;

    VEC{k,1}.U =  -A*pi*sin(pi*f).*cos(pi*VEC{k}.Y);
    VEC{k,1}.V = A*pi*(2*a*VEC{k}.X + b).*cos(pi*f).*sin(pi*VEC{k}.Y);
    VEC{k,1}.C = logical(ones(length(y),length(x)));
end

%% Animate Flow
if nargin == 3  && strcmpi(varargin{1},'Animate')
    FlowType = 'Double Gyre';
    sky = 2;
    %     skx = round((20/6)*sky);
    skx = 2;
    %     h = figure('Units','normalized','Position',[0 0 3.3*Height Height])
%     figure('visible','on')
    for k = 1:5:length(VEC)
        quiver(VEC{1}.X(1:sky:end,1:skx:end),VEC{1}.Y(1:sky:end,1:skx:end),...
            VEC{k}.U(1:sky:end,1:skx:end),VEC{k}.V(1:sky:end,1:skx:end),'k')
        axis tight
        title([FlowType,',Time=',num2str(round(t(k),2))],'interpreter','latex','FontSize',14);
        pbaspect([2 1 1])
        pause(1/length(VEC));
        %         drawnow limitrate
    end
    %     close(h)
end