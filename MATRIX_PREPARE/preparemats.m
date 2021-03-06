clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/FEM/')

conf = '3';

for conf={'1', '2', '3', '4', '4a', 'D1'}
  conf = conf{1};
fname = sprintf(['../../ABAQUS_MODELS_TWOS/CONF%s/PROCESS/JOBRUN/' ...
                 'Modelmats.mat'], conf);
load(fname, 'M', 'K', 'R', 'Fv');
Fv = Fv';
R = R';

%% Extract Mesh
MESH.Nds  = dlmread(sprintf('../../ABAQUS_MODELS_TWOS/CONF%s/PROCESS/Nodes.dat', conf));
MESH.Quad = dlmread(sprintf('../../ABAQUS_MODELS_TWOS/CONF%s/PROCESS/Elements.dat', conf));
MESH.Quad(:, 2:end) = MESH.Quad(:, 2:end)+1;
MESH.Tri = [];

MESH.Nn = size(MESH.Nds,1);
MESH.Ne_Quad = size(MESH.Quad,1);
MESH.Ne_Tri = size(MESH.Tri,1);
MESH.Ne = MESH.Ne_Quad + MESH.Ne_Tri;
MESH.dpn = 3;

[Q1, T1] = ZTE_ND2QP(MESH, 1);
%% Patch Selection
xmin_max_rng = [min(MESH.Nds(:,1)) max(MESH.Nds(:,1)) range(MESH.Nds(:,1))];
ymin_max_rng = [min(MESH.Nds(:,2)) max(MESH.Nds(:,2)) range(MESH.Nds(:,2))];

E_ctrds = [mean(reshape(MESH.Nds(MESH.Quad(:, 2:end), 1), [], 4), 2), ...
	   mean(reshape(MESH.Nds(MESH.Quad(:, 2:end), 2), [], 4), 2)];  % Element Centroids

if ismember(conf, {'1', '2', '3', '4', '4a'})  # Square interface
  Npatches = 4;  % Must be a square number
  
  xs = linspace(xmin_max_rng(1), xmin_max_rng(2), sqrt(Npatches)+1);
  ys = linspace(ymin_max_rng(1), ymin_max_rng(2), sqrt(Npatches)+1);

  Area = cell(Npatches, 1);
  qps = cell(Npatches, 1);
  Pnds = cell(Npatches, 1);
  Pels = cell(Npatches, 1);
  Ctrds = zeros(Npatches, 2);
  for i=1:sqrt(Npatches)
    for j=1:sqrt(Npatches)
      p = (i-1)*sqrt(Npatches) + j;

      Pels{p} = find(E_ctrds(:,1) >= xs(i) & E_ctrds(:,1) <= xs(i+1) & ...
		     E_ctrds(:,2) >= ys(j) & E_ctrds(:,2) <= ys(j+1));
      Pnds{p} = unique(reshape(MESH.Quad(Pels{p}, 2:end), [], 1));
      Area{p} = full(sum(T1(Pnds{p}, Pels{p}), 2));
      qps{p} = Q1(Pels{p}, :)*MESH.Nds;
      Ctrds(p, :) = sum(repmat(Area{p}, 1,2).*MESH.Nds(Pnds{p}, 1:2))/sum(Area{p});
    end
  end
else
  Npatches = 5;
  xs = linspace(xmin_max_rng(1), xmin_max_rng(2), Npatches+1);

  Area = cell(Npatches, 1);
  qps = cell(Npatches, 1);
  Pnds = cell(Npatches, 1);
  Pels = cell(Npatches, 1);
  Ctrds = zeros(Npatches, 2);
  for p=1:Npatches
      Pels{p} = find(E_ctrds(:,1) >= xs(p) & E_ctrds(:,1) <= xs(p+1));
      Pnds{p} = unique(reshape(MESH.Quad(Pels{p}, 2:end), [], 1));
      Area{p} = full(sum(T1(Pnds{p}, Pels{p}), 2));
      qps{p} = Q1(Pels{p}, :)*MESH.Nds;
      Ctrds(p, :) = sum(repmat(Area{p}, 1,2).*MESH.Nds(Pnds{p}, 1:2))/sum(Area{p});
  end
end
PatchAreas = full(cellfun(@(c) sum(c), Area));

%% Relative Coordinate Transformation: [XT-XB; XB; eta]
Ngen = size(K,1)-2*MESH.Nn*MESH.dpn
Trel = [kron([1 1; 0 1], speye(MESH.Nn*MESH.dpn)), sparse(2*MESH.Nn*MESH.dpn, Ngen);
	zeros(Ngen, 2*MESH.Nn*MESH.dpn) speye(Ngen)];
Krel  = Trel'*K*Trel;  Krel = 0.5*(Krel+Krel');
Mrel  = Trel'*M*Trel;  Mrel = 0.5*(Mrel+Mrel');
Fvrel = Trel'*Fv;
Rrel  = R*Trel;

%% Weak Form Integral Matrices
NTN = sparse(MESH.Nn*3, MESH.Nn*3);
NTG = sparse(MESH.Nn*3, Npatches*6);
GTG = sparse(Npatches*6, Npatches*6);
BNV = sparse(MESH.Nn, Npatches);  % Force Transformation Matrix
for n=1:Npatches
    [P, Nums, NTNmat, NTGmat, GTGmat] = CONSPATCHMAT(MESH.Nds, [], MESH.Quad(Pels{n}, :), Ctrds(n, :));
    
    NTN = NTN + NTNmat;
    NTG(:, (n-1)*6+(1:6)) = NTG(:, (n-1)*6+(1:6)) + NTGmat;
    GTG((n-1)*6+(1:6), (n-1)*6+(1:6)) = GTG((n-1)*6+(1:6), (n-1)*6+(1:6)) + GTGmat;
    
    BNV(Pnds{n}, n) = 1.0/sum(Area{n});
end
GN = GTG\(NTG');

%% Reduced Model
ni = 1:(MESH.Nn*MESH.dpn);
mi = setdiff(1:length(Krel), ni);
cnum = 1e10;

Ka = sparse([Krel(ni, ni)+cnum*NTN Krel(ni, mi) -cnum*NTG;
             Krel(mi, ni), Krel(mi, mi), zeros(length(mi),Npatches*6);
             -cnum*NTG' zeros(Npatches*6, length(mi)) cnum*GTG]);
Ka = 0.5*(Ka+Ka');
Ma = sparse(blkdiag(Mrel, zeros(Npatches*6)));
Ma = 0.5*(Ma+Ma');
Fva = [Fvrel; zeros(Npatches*6,1)];
LamTa = [NTG*inv(GTG); zeros(length(mi)+Npatches*6,Npatches*6)];
Ftmpa = Ka*sparse([zeros(length(Krel), 1); cnum^(-1)*NTG\Fvrel(ni)]);

Ra = sparse([Rrel zeros(size(Rrel,1), Npatches*6)]);

dofred = 6;
vdofs = length(K) + reshape(1:(Npatches*6), 6, Npatches);
vdofs = reshape(vdofs(1:dofred,:), Npatches*dofred, 1);  % virtual DOFs

[Mhcb, Khcb, Thcb] = HCBREDUCE(Ma, Ka, vdofs, Ngen);
Mhcb = 0.5*(Mhcb+Mhcb');  Khcb = 0.5*(Khcb+Khcb');
Fvhcb = Thcb'*Fva;
Rhcb = Ra*Thcb;
Lamhcb = Thcb'*LamTa;
disp('DONE!')

%% NULL-SPACE REDUCTION
[Vhcb, Dhcb] = eig(full(Khcb(Npatches*dofred+1:end, Npatches*dofred+1:end)),
		   full(Mhcb(Npatches*dofred+1:end, Npatches*dofred+1:end)));
[Dhcb, si] = sort(diag(Dhcb));  Vhcb = Vhcb(:, si);
Vrbms = [zeros(Npatches*dofred, 6); Vhcb(:, 1:6)];  Vrbms = Vrbms./sqrt(diag(Vrbms'*Mhcb*Vrbms)');
L = null(Vrbms'*Mhcb);

M = L'*Mhcb*L;   M = 0.5*(M+M');
K = L'*Khcb*L;   K = 0.5*(K+K');
Fv = L'*Fvhcb;
R = Rhcb*L;
Th = Thcb*L;
LamT = L'*Lamhcb;
[sort(eig(K), 'descend') sort(eig(M), 'descend')]
disp(cnum)

%% %% Check Eigenvalues
%% [V, D] = eigs(K, M, 30, 'SM');
%% [Vrel, Drel] = eigs(Krel, Mrel, 30, 'SM');

%% Save
fname = sprintf('MAT_NULLRED_C%s.mat', conf);
save(fname, 'M', 'K', 'L', 'Fv', 'R', 'Th', 'LamT', 'MESH', ...
     'PatchAreas', 'Npatches', 'dofred', 'NTN', 'GTG', 'NTG', ...
     'cnum', 'Pels', 'Pnds')
disp('SAVED!');
end
