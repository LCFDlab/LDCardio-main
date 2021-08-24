
% By Ahmed Darwish
% Last update : 11 May 2021

function [varargout] = LG_Descriptors (A_FWD,A_BWD,P,dt,varargin)

%% Reshape the A_FWD and B_WD into matrix form 
% The "if" specifies either the 'VIN-DLD' ( variable iteration number DLD)
% algorithm is used in the computations or the regular algorithm (constant
% number of iterations). See Ref [1] for more details about the VIN-DLD
% algorithm.

if nargout>5 && nargin==3
    error('Number of output arguments exceeds 5.\n Please specify VIN-DLD',...
        'or revise the list of output')
end

if nargin == 5 && strcmpi(varargin{1},'VIN-DLD')
    for k = 1 : size(A_FWD.X,2)
        
        if (any(isnan(A_BWD.X(:,k))))||(any(isnan(A_FWD.X(:,k))))
            if (any(isnan(A_BWD.X(:,k))))
                Indx_bwd(k) = find(isnan(A_BWD.X(:,k)), 1, 'first');
            else
                Indx_bwd(k) = size(A_BWD.X,1);
            end
            if (any(isnan(A_FWD.X(:,k))))
                Indx_fwd(k) = find(isnan(A_FWD.X(:,k)), 1, 'first');
                
            else
                Indx_fwd(k) = size(A_FWD.X,1);
            end
            Indx(k) = min([Indx_fwd(k) Indx_bwd(k)]);
            A_BWD.X(Indx(k)+1:end,k) = zeros(size(A_BWD.X,1)-Indx(k),1);
            A_BWD.Y(Indx(k)+1:end,k) = zeros(size(A_BWD.X,1)-Indx(k),1);
            A_FWD.X(Indx(k)+1:end,k) = zeros(size(A_BWD.X,1)-Indx(k),1);
            A_FWD.Y(Indx(k)+1:end,k) = zeros(size(A_BWD.X,1)-Indx(k),1);
        else
            Indx_fwd(k) = size(A_FWD.X,1);
            Indx_bwd(k) = size(A_BWD.X,1);
        end
    end
    
    EuclX_FWD = abs(A_FWD.X(2:end,:)-A_FWD.X(1:end-1,:)).^P;
    EuclY_FWD = abs(A_FWD.Y(2:end,:)-A_FWD.Y(1:end-1,:)).^P;
    EuclX_BWD = abs(A_BWD.X(2:end,:)-A_BWD.X(1:end-1,:)).^P;
    EuclY_BWD = abs(A_BWD.Y(2:end,:)-A_BWD.Y(1:end-1,:)).^P;
    
    FWD_DLD = sum(((EuclX_FWD + EuclY_FWD)).^(1/P),1,'omitnan')./(Indx_fwd(k)*dt);
    BWD_DLD = sum(((EuclX_BWD + EuclY_BWD)).^(1/P),1,'omitnan')./(Indx_bwd(k)*dt);
    
    %     for k = 1 : size(A_FWD.X,2)
    %
    %         FWD_DLD(k) = sum(((EuclX_FWD(:,k) + EuclY_FWD(:,k))).^(1/P),1,'omitnan');
    %         BWD_DLD(k) = sum(((EuclX_BWD(:,k) + EuclY_BWD(:,k))).^(1/P),1,'omitnan');
    %
    %     end
    DLD = FWD_DLD + BWD_DLD;
    % List of possible output arguments [NOTICE THEIR ORDER !!!]
    varargout{1} = BWD_DLD;
    varargout{2} = FWD_DLD;
    varargout{3} = DLD;
    varargout{4} = A_BWD;
    varargout{5} = A_FWD;
    varargout{6} = Indx_bwd;
    varargout{7} = Indx_fwd;
    
    
    % ----- If the regular DLD algorithm is selected
else
    
    EuclX_FWD = abs(A_FWD.X(2:end,:)-A_FWD.X(1:end-1,:)).^P;
    EuclY_FWD = abs(A_FWD.Y(2:end,:)-A_FWD.Y(1:end-1,:)).^P;
    EuclX_BWD = abs(A_BWD.X(2:end,:)-A_BWD.X(1:end-1,:)).^P;
    EuclY_BWD = abs(A_BWD.Y(2:end,:)-A_BWD.Y(1:end-1,:)).^P;
    
    FWD_DLD = sum(((EuclX_FWD + EuclY_FWD)).^(1/P),1,'omitnan');
    BWD_DLD = sum(((EuclX_BWD + EuclY_BWD)).^(1/P),1,'omitnan');
    
    %     for k = 1 : size(A_FWD.X,2)
    %         FWD_DLD(k) = sum(((EuclX_FWD(:,k) + EuclY_FWD(:,k))).^(1/P),1,'omitnan');
    %         BWD_DLD(k) = sum(((EuclX_BWD(:,k) + EuclY_BWD(:,k))).^(1/P),1,'omitnan');
    %     end
end
DLD = FWD_DLD + BWD_DLD;

% List of possible output arguments [NOTICE THEIR ORDER !!!]
varargout{1} = BWD_DLD;
varargout{2} = FWD_DLD;
varargout{3} = DLD;
varargout{4} = A_BWD;
varargout{5} = A_FWD;
end
% -----NOTES -------
% In lines 56-60 and 90-94 the 2-norm can be computed using the vecnorm
% however the use of vecnorm will consider NaN values. Therefore, the code
% computes the 2-norm with sum and 'omitnan' to remove the NaN elements.

% ----REFERENCES---------
% [1] - Garcıa-Garrido, 2020, An Extension of Discrete Lagrangian Descriptors
% for Unbounded Maps, International Journal of Bifurcation and Chaos.