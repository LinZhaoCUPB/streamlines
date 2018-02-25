function [sl, TOF] = streamline_pollock(Ic, If, Ixyz, G, vt, proc, poro)
% Traces individual streamline for the cartGrid based on the Pollock's (1988) algorithm
%
% PARAMETERS:
% Ic    -  start cell
% If    -  start face
% Ixyz  -  start coordinates
% G     -  grid structure
% vt    -  face velocity
% proc  -  producer cells
% poro  -  porosity
%
% RETURNS:
% sl    -  coordinates of single streamline
% TOF   -  time of flight of each cell
%
% Written by Lin Zhao, CUPB, CHINA
%
% SEE ALSO:
%   `generate_streamlines`, `streamline_ck`


% ----------------------
N  = double(G.faces.neighbors);
DXYZ = G.DXYZ; % spatial step for CartGrid
% -------------------------
sl      = zeros(1,3);
TOF     = zeros(1,2);
sl(1,:) = Ixyz;
% main loop ---------------------------------
k = 2;
while true
    Ec = N(If,:);
    Ec = Ec(Ec~=Ic); % exit cell

    if any(Ec == proc);break;end % reaches producer cells, generation done
    
    Ef = getfaces(G,Ec);    % possible exit faces
    
    vEf = vt(Ef);    % face velocity
    
    vEf = reshape(vEf, 2, 3)'; % reshape by [x- x+; y- x+; z- z+]

    dxyz = DXYZ{Ec}'; % [dx dy dz]'
    dv = (vEf(:,2) - vEf(:,1))./dxyz; % velocity gradient
    
    move = G.cells.centroids(Ec,:) - dxyz'/2; % move cell to origin
    
    Ixyz0 = Ixyz - move;
    
    vI = dv .* Ixyz0' + vEf(:,1); % velocity of inlet point
    
    vI2 = [vI, vI]; % extend for computing tof
    dv2 = [dv, dv]; % extend for computing tof
    
    tofs = log(vEf./vI2)./dv2; % possible tofs
    
    tofs = reshape(tofs', [], 1);
    tofs(tofs <= 0  )  = inf;
    tofs(Ef == If   )  = inf;
    tofs(isnan(tofs))  = inf;
    
    tof = min(tofs); % min tof is selected 
    
    Ef  = Ef(tofs == tof); % exit face
    
    Exyz0 = (exp(dv.*tof).*vI - vEf(:,1))./dv; % exit coordinate0
    
    % modify exit coordinate based on face coordinate
    box = [zeros(3,1), dxyz];
    box = reshape(box', [], 1);
    minind = find(tofs == tof);
    Exyz0(ceil(minind/2)) = box(minind); 
    
    Exyz  = Exyz0' + move; % exit coordinate 
    
    if any(Exyz0 < 0); 
        error('exit coordinate < 0')
    end
    
    % output ---------------------------------
    tof = tof.*poro(Ec); % real tof
    TOF(k-1,1) = tof;
    TOF(k-1,2) = Ec;
    sl(k,:)    = Exyz;
    
    % next cell ---------------------------------
    If   = Ef;
    Ic   = Ec;
    Ixyz = Exyz;
    k = k + 1;
end
end

%%
function faces = getfaces(G, cell)
% sorted faces of cell
% face 1 -- x -
% face 2 -- x +
% face 3 -- y -
% face 4 -- y +
% face 5 -- z -
% face 6 -- z +

facePos = G.cells.facePos(cell) : G.cells.facePos(cell+1)-1;

faces = G.cells.faces(facePos, 1);
ind   = G.cells.faces(facePos, 2);

[~, i] = sort(ind);
faces  = faces(i);
end
