conf = 'D1';
fname = sprintf('MAT_NULLRED_%s.mat', conf);

load(fname);
save(fname, 'M', 'K', 'L', 'Fv', 'R', 'Th', 'LamT', 'MESH', ...
     'PatchAreas', 'Npatches', 'dofred', 'NTN', 'GTG', 'NTG', ...
     'cnum', 'Pels', 'Pnds', '-7')
