function [sc, sl, tof] = generate_streamlines(G, W, qt, q0, poro)
% Generates streamlines for a flow field
%
% PARAMETERS:
% G     grid structure, Cart or Corner-Point grid
% W     well structure
% qt    face flux
% q0    flux on each streamline
% poro  porosity
%
% RETURNS:
% sc    start configuration
% sl    coordinates of all streamlines
% TOF   time of flight of all streamlines
%
% Written by Lin Zhao, CUPB, CHINA
%
% SEE ALSO:
%   `streamline_pollock`, `streamline_ck`

% ----------------------
% well cells and directions
wc    = vertcat(W.cells);
dir   = vertcat(W.dir);

sign = vertcat(W.sign);
sign = rldecode(sign, cellfun(@numel, {W.cells}));

cstatus  = vertcat(W.cstatus);

wc     = wc  (cstatus);
dir    = dir (cstatus);
sign   = sign(cstatus);

% x - char(120) y - char(121) z - char(122)
dir = double(dir) - 119;

proc   = wc (sign == -1);  % producer cells
injc   = wc (sign ==  1);  % injector cells
injdir = dir(sign ==  1);  % injector directions

% ----------------------
% start configuration
% 1 - start cell; 2 - start face; 3:5 - start coordinates
sc = arrayfun(@(cell, dir)start_configuration(G, qt, q0, cell, dir),...
    injc, injdir,  'UniformOutput', false);
sc = cell2mat(sc);
%%
if strcmp(G.type{end}, 'corner-point grid')
    fun_streamline = @streamline_ck;
    % face flux
    xt = qt;
elseif strcmp(G.type{end}, 'cartesian grid')
    fun_streamline = @streamline_pollock;
    % face velocity
    xt = qt./G.faces.areas;
    % spatial step for CartGrid
    G.DXYZ = arrayfun(@(cell)get_spatialstep(G, cell), (1:G.cells.num)', 'UniformOutput', false);
else
end
%% streamline generation
sln = size(sc,1);
[sl, tof]  = deal(cell(sln,1));

h = waitbar(0, sprintf('Generating %4d of %4d Streamlines', 0 , sln));

for s = 1:sln
    waitbar(s/sln, h, sprintf('Generating %4d of %4d Streamlines', s , sln));
    try
        [sl{s,1}, tof{s,1}] = ...
            fun_streamline(sc(s,1), sc(s,2), sc(s,3:5), G, xt, proc, poro);
    catch
    end
end

waitbar(1, h, 'Generation Done!')
%%
% remove error streamlines
err = cellfun(@isempty,sl);
sc(err,:) = [];
sl(err)   = [];
tof(err)  = [];
end

%%
function sc = start_configuration(G, qt, q0, cell, dir)
% start configuration
% 1 - start cell; 2 - start face; 3:5 - start coordinates
faces = getfaces(G, cell);

% start faces, based on well direction
faces  = faces(~ismember(1:6, [dir*2-1 dir*2]));

% middile points of two pillars
pmc = arrayfun(@(ii)pillar_mid_coords(G, faces, ii), 1:4, 'UniformOutput', false);

% streamline number
sl_num = abs(qt(faces)/q0);
sl_num = ceil(sl_num);

% start point gap
fgap    = @(x,y,num)(y-x)/(num+1);
fpoints = @(x,y,num)bsxfun(@plus, x, (1:num)' * fgap(x, y, num));

% start points
points = cellfun(@(p, num)fpoints(p(1,:), p(2,:), num), pmc', num2cell(sl_num), 'UniformOutput', false);
points = cell2mat(points);

% extend faces for storage
faces_all = arrayfun(@(x,num)x*ones(num,1), faces, sl_num, 'UniformOutput', false);
faces_all = cell2mat(faces_all);

% extend cell for storage
cell_all  = cell*ones(size(faces_all));

% start configuration
sc = [cell_all, faces_all, points];
end

%%
function pmc = pillar_mid_coords(G, faces, ii)
fnodePos = @(f)G.faces.nodePos(f) : G.faces.nodePos(f+1)-1;
fnodes   = @(f)G.faces.nodes(fnodePos(f));

% selected face
face_i    = faces(ii);
% rest faces
faces_r   = faces; faces_r(ii) = [];

nodes_i   =  fnodes(face_i);
nodes_r    = arrayfun(fnodes, faces_r, 'UniformOutput', false);

% pillar is computed by face intersection
pillar_nodes   = cellfun(@(x)intersect(nodes_i, x), nodes_r, 'UniformOutput', false);
pillar_nodes   = pillar_nodes(~cellfun(@isempty, pillar_nodes));

% middile coordinates of two pillars
pr_mid_coords  = cellfun(@(x)mean(G.nodes.coords(x,:),1), pillar_nodes, 'UniformOutput', false);
pmc = cell2mat(pr_mid_coords);
end

%%
function ss = get_spatialstep(G, cell)
% spatialstep (dx dy dz) for each grid
faces = getfaces(G, cell);
fc    = G.faces.centroids(faces,:);
ss    = cellfun(@(x)max(x)-min(x), num2cell(fc,1));
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
%%

