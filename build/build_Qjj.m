function Qjj = build_Qjj(Ma,k_red,c_ref,Pa,SYM)
% build_Qjj   Aerodynamic Influence Coeffficient AIC Matrix 
% INPUT
%   Ma       : Mach number
%   k_red    : Vector of reduced frequencies k_red = omega*c_ref/(2*V_inf)
%   c_ref    : Reference Chord for normalization
%   Pa       : cell array, one entry per aerodynamic panel in aerodynamic
%              coordinate system with freestream in positive x-dir.
%              For each panel the four corner nodes are stored as
% 
%               P{panel}{node} = [ x_pos;         1---4
%                                  y_pos;         |   |
%                                  z_pos ]        2---3
% 
% SYM        : Symmetrie Flag Wing mirrored at xz-plane
%
% OUTPUT
% Qjj        : complex valued AIC Matrix [nj,nj,nk] from dimensionless downwash to panel differential pressure coefficient
%              depends on Ma, k_red = omega * c_ref / ( V_inf * 2 ) and the panel geometry P
%              maps sinusoidally oscillating downwash to sinusiodally oscillating defferential pressure coefficients (amplitude und phase)
%              Qjj=inv(Ajj)
%              cpj = Ajj^-1 * wj
%
%              wj = u_zj/V_inf        Normalized (by freestream vel) downwash at pistolesi points     
%              cpj = Delta_pj/q_bar   Normalized (by dynamic pressure) pressure difference at the j-th panel
%
% Sign convention
% _________________________________________________________________________
% Freestream direction is positiv x-axis in aero coord. sys.
%
% Jonas * July 2025
% _________________________________________________________________________


nj = length(Pa);
nk = length(k_red);

if SYM == 1
    Pa_sym = cell(2*nj,1);
    for i_p = 1:nj
        for i_k = 1:4
            Pa_sym{i_p}{i_k}= Pa{i_p}{i_k};
        end
        Pa_sym{i_p+nj}{4} = Pa_sym{i_p}{1};
        Pa_sym{i_p+nj}{3} = Pa_sym{i_p}{2};
        Pa_sym{i_p+nj}{2} = Pa_sym{i_p}{3};
        Pa_sym{i_p+nj}{1} = Pa_sym{i_p}{4};
        for i_k = 1:4
            Pa_sym{i_p+nj}{i_k}(2) = -Pa_sym{i_p+nj}{i_k}(2);
        end
    end
    Pa=Pa_sym;
    np=2*nj;
else
    np=nj;
end

BVP1     = zeros(3, np);
BVP2     = zeros(3, np);
JP       = zeros(3, np);
LP       = zeros(3, np);
N        = zeros(3, np);
c_mean   = zeros(1, np);
s        = zeros(1, np);
for i_p = 1:np
    c1 = abs(Pa{i_p}{1}(1) - Pa{i_p}{2}(1));
    c2 = abs(Pa{i_p}{4}(1) - Pa{i_p}{3}(1));
    c_mean(i_p)     = 0.5*(c1+c2);

    JP1 = 0.75*Pa{i_p}{2} + 0.25*Pa{i_p}{1};
    JP2 = 0.75*Pa{i_p}{3} + 0.25*Pa{i_p}{4};
    JP(:,i_p)  = 0.5*(JP1+JP2);

    BVP1(:,i_p)     = 0.25*Pa{i_p}{2} + 0.75*Pa{i_p}{1};
    BVP2(:,i_p)     = 0.25*Pa{i_p}{3} + 0.75*Pa{i_p}{4};
    LP(:,i_p) = 0.5*(BVP1(:,i_p)+BVP2(:,i_p));

    % s(i_p) = norm(Pa{i_p}{4}(2:3) - Pa{i_p}{1}(2:3));
    S = Pa{i_p}{4}-Pa{i_p}{1} - (Pa{i_p}{2}-Pa{i_p}{1})/norm((Pa{i_p}{2}-Pa{i_p}{1}))^2*dot((Pa{i_p}{2}-Pa{i_p}{1}),(Pa{i_p}{4}-Pa{i_p}{1}));
    s(i_p) = norm(S);
    % S = S/s;
    
    Np=cross((Pa{i_p}{2} - Pa{i_p}{1}), (Pa{i_p}{4}-Pa{i_p}{1}));
    Np=Np/norm(Np);
    N(:,i_p)=Np;
end

pa.CP = JP';           % (:,3)     % Control point on 3/4 line
pa.SP = LP';           % (:,3)     % Sending point on 1/4 line midspan
pa.D1 = BVP1';         % (:,3)     % inner border point on doublet line (1/4 line)
pa.D5 = BVP2';         % (:,3)     % outer border point on doublet line
pa.N =  N';            % (:,3)     % normal vector
pa.chord = c_mean';    % double    % average box chord
pa.s = s';             % double    % span of box

int      = "Quartic";   % "Parabolic"; 
app      = "Desmarais"; % "Watkins"; "Laschka";
pkern    = 0;           % plot Kernel
dsp      = 0;           % verbose


Qjj = zeros(nj,nj,nk);
for i_k=1:nk
    kr = k_red(i_k);
    kred = kr*2/c_ref; % specific reduced frequency used in DLM [ omega/V_inf ]
    D_vlm = VLM(pa,Ma,dsp); % Get VLM (steady state effect)
    D_dlm = DLM(pa,Ma,kred,int,app,pkern,dsp); % Get DLM (oscillatory effect)

    % Combine VLM and DLM matrices
    D = -(D_vlm+D_dlm);

    % Concatenate results to one wing
    if SYM == 1
        D_rr = D(1:nj,1:nj);        % influence from right wing on itself
        D_rl = D(1:nj,nj+1:2*nj);   % influence from left wing on right wing
        D = D_rr+D_rl;
    end
    Qjj(:,:,i_k) = inv(D);          % AIC Matrix DCp_j = Qjj * w_j
end

end