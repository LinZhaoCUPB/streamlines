function [sl, TOF] = streamline_ck(Ic, If, Ixyz, G, qt, proc, poro, varargin)
% Traces individual streamline for the corner-point Grid based on the Cordes and Kinzerlbach's (1992) algorithm
%
% PARAMETERS:
% Ic    -  start cell
% If    -  start face
% Ixyz  -  start coordinates
% G     -  grid structure
% qt    -  face flux
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
%   `generate_streamlines`, `streamline_pollock`

% ---------------------
opt = struct('real_tof', false);
opt = merge_options(opt, varargin{:});
% ---------------------
N  = double(G.faces.neighbors);
% -------------------------
sl      = zeros(0,3);
TOF     = zeros(0,2);
sl(1,:) = Ixyz;
% main loop ---------------------------------
k = 2;
while true
    Ec = N(If,:);    
    Ec = Ec(Ec~=Ic); % exit cell
    
    if any(Ec == proc);break;end % reaches producer cells, generation done
    
    Ef = getfaces(G, Ec);    % possible exit faces
    
    Iabc  = trilinear_interp_inverse(G, Ixyz, Ec);  % inlet abc coordinates (unit system)
    
    Iabc(Iabc<0) = 0; Iabc(Iabc>1) = 1;  % modify abc coordinates
    
    qEf = qt(Ef); % face flux
    
    qEf = reshape(qEf, 2, 3)'; % reshape by [x- x+; y- x+; z- z+]
    
    dq = qEf(:,2) - qEf(:,1);  % flux gradient
    m  = qEf(:,1);
    
    qI = dq .* Iabc' + m;  % flux of inlet point
    
    qI2 = [qI, qI]; % extend for computing tof
    dq2 = [dq, dq]; % extend for computing tof
    
    tofs = log(qEf./qI2)./dq2; % possible tofs
    
    tofs = reshape(tofs', [], 1);
    tofs(tofs <= 0  )  = inf;
    tofs(Ef == If   )  = inf;
    tofs(isnan(tofs))  = inf;
    
    tof = min(tofs); % min tof is selected
    
    Ef  = Ef(tofs == tof); % exit face
    
    Eabc = (exp(dq.*tof).*qI - m)./dq; % exit abc coordinate
    
    % modify exit coordinate based on face coordinate
    box = [zeros(3,1), ones(3,1)];
    box = reshape(box', [], 1);
    minind = find(tofs == tof);
    Eabc(ceil(minind/2)) = box(minind); 
    
    Exyz  = trilinear_interp(G, Eabc, Ec); % exit xyz coordinate
    
    % output ---------------------------------
    
    % compute tof in xyz coordinates by integration,but it may take some
    % time ...
    if opt.real_tof
        tof  = tof_by_integJac(G, Ec, qI, dq, m, tof);
    end
    
    tof  = tof.*poro(Ec);
    TOF(k-1,1) = tof;
    TOF(k-1,2) = Ec;
    sl(k,:)    = Exyz;
    
    % next cell ---------------------------------
    If   = Ef;
    Ic   = Ec;
    Ixyz = Exyz;
    k    = k + 1; 
end
end

%%
function Ixyz  = trilinear_interp(G, Iabc, cell)
% based on Tri-linear Interpolation
% slove the  x y z coordinate from a b c coordinate
coords = get_stand_coords(G, cell);
N = initN;
P = N*coords;
a = Iabc(1); b = Iabc(2); c = Iabc(3); 
s = [a, b, c, a*b, b*c, a*c, a*b*c, 1] * P;
Ixyz = s;
end

%%
function Iabc  = trilinear_interp_inverse(G, Ixyz, cell)
% based on Tri-linear Interpolation
% slove the  a b c coordinate from x y z coordinate
coords = get_stand_coords(G, cell);
N = initN;
P = N*coords;
% since the expression of a b c is implicit, we solve them by AD approach
% initialization
a = 0; b = 0; c = 0;
dx = [inf inf inf];
tol = 1e-4;
while max(abs(dx)) > tol
    [a,b,c] = initVariablesADI(a,b,c);
    eqs = P' * [a; b; c; a*b; b*c; a*c; a*b*c; 1] - Ixyz';
    J =  -cell2mat(eqs.jac);
    v =  eqs.val;
    dx = J\v;
    a = double(a) + dx(1);
    b = double(b) + dx(2);
    c = double(c) + dx(3);
end
Iabc = [a, b, c];
end

%%
function tof  = tof_by_integJac(G, cell, qI, dq, bq, tof)
% compute tof by integration
coords = get_stand_coords(G, cell);
N = initN;
P = N*coords;

syms a b c qia qib qic ca cb cc ba bb bc t;
s = [a, b, c, a*b, b*c, a*c, a*b*c, 1] * P;

% jacobian matrix
jacdet = det(jacobian(s,[a,b,c]));
f_jacdet = @(a,b,c)eval(jacdet);

at = (qia*exp(ca*t)-ba)/ca;
bt = (qib*exp(cb*t)-bb)/cb;
ct = (qic*exp(cc*t)-bc)/cc;

% interpolation from 0 to t in abc coordinate
tof_syms= int(f_jacdet(at, bt, ct), t, 0, t);

f_tof  = @(qia, qib, qic, ca, cb, cc, ba, bb, bc, t)eval(tof_syms);

% tof
tof = f_tof(qI(1), qI(2),qI(3),dq(1),dq(2),dq(3),bq(1),bq(2),bq(3),tof);

if tof <= 0
    error('tof <= 0')
end
end
%%
function N = initN
% mapping matrix
N =zeros(8,8);
N(1,2) =  1;  N(1,1) = -1;
N(2,4) =  1;  N(2,1) = -1;
N(3,5) =  1;  N(3,1) = -1;
N(4,1) =  1;  N(4,3) =  1;  N(4,2) = -1;  N(4,4) = -1;
N(5,1) =  1;  N(5,8) =  1;  N(5,4) = -1;  N(5,5) = -1;
N(6,1) =  1;  N(6,6) =  1;  N(6,2) = -1;  N(6,5) = -1;
N(7,2) =  1;  N(7,4) =  1;  N(7,5) =  1;  N(7,7) =  1;
N(7,1) = -1;  N(7,3) = -1;  N(7,6) = -1;  N(7,8) = -1;
N(8,1) =  1;
end
%%
function scoords = get_stand_coords(G, cell)
% get standard coordinates, used for corner-point grid
% 1 -- x0 y0 z0
% 2 -- x1 y0 z0
% 3 -- x1 y1 z0
% 4 -- x0 y1 z0
% 5 -- x0 y0 z0
% 6 -- x0 y0 z0
% 7 -- x0 y0 z0
% 8 -- x0 y0 z0
facePos = getfaces(G,cell);
f_nodePos = @(f)G.faces.nodes(G.faces.nodePos(f):G.faces.nodePos(f+1)-1,:);
f_node = @(k)G.nodes.coords(f_nodePos(facePos(k)),:);
% get coordinates from face intersections
% face 1 -- x -
% face 2 -- x +
% face 3 -- y -
% face 4 -- y +
% face 5 -- z -
% face 6 -- z +
% node 1 -- intersect face 1 3 5
% node 2 -- intersect face 2 3 5
% node 3 -- intersect face 2 4 5
% node 4 -- intersect face 1 4 5
% node 5 -- intersect face 1 3 6
% node 6 -- intersect face 2 3 6
% node 7 -- intersect face 2 4 6
% node 8 -- intersect face 1 4 6
f_intersect = @(x)intersect(f_node(x(1)), intersect(f_node(x(2)),f_node(x(3)),'rows'), 'rows');
T =[1 3 5;2 3 5;2 4 5;1 4 5;1 3 6;2 3 6;2 4 6;1 4 6];

scoords = cellfun(f_intersect, num2cell(T,2), 'UniformOutput', false);
scoords = cell2mat(scoords);

% for j = 1:8
%     coord(j,:) = f_intersect(T(j,:));
% end
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