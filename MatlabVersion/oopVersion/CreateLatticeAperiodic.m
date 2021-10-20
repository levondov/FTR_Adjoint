% define a lattice
function obj = CreateLatticeAperiodic(obj, a, period, numPeriodicCells, runMoments)

if nargin < 5
    runMoments = false;
end

drift = 0.04;
db = 0.15;
qlength = 0.04;
qlengthSpecial = 0.04;
dbSpecial = 0.15 - db*0.05; % Pertrubation

% special starting quads
qstart1 = [0.0, qlengthSpecial/2.0+drift, qlengthSpecial*1.5+drift*2];
db1 = [dbSpecial, -dbSpecial, db];
qend1 = [qstart1(1)+qlengthSpecial/2.0, qstart1(2)+qlengthSpecial, qstart1(3)+qlength/2.0];
qrot1 = [0.0,0.0,0.0];

% regular lattice
qstart = [0.0, 0.06, 0.14];
db = [0.15*a(1), -0.15*a(3), 0.15];
qend = [qstart(1)+qlength/2.0, qstart(2)+qlength, qstart(3)+qlength/2.0];
qrot = [0.0, 0.0, 0.0];

% create periodic lattice
lattice1 = [qstart1',qend1',db1',qrot1'];
lattice2 = [qstart',qend',db',qrot'];
lattice2(:,[1,2]) = lattice2(:,[1,2]) + qend1(3);

% repeating lattice elements
latticeCopy = lattice2;
zend = (qend(3)-qstart(1));
for i = 1:(numPeriodicCells-1)
    latticeCopy(:,[1,2]) = latticeCopy(:,[1,2]) + zend;
    lattice2 = [lattice2; latticeCopy];
end

lattice = [lattice1;lattice2];

% update last element to match first one
lattice(end,2) = lattice(end,1) + qlengthSpecial/2.0;
lattice(end,3) = dbSpecial;
lattice(end,4) = qrot1(1);

db = lattice(:,3);
qstart = lattice(:,1);
qend = lattice(:,2);
qrot = lattice(:,4);
zend = lattice(end,2);
zstart = lattice(1,1);

obj = obj.CreateLatticeProfile(db,qstart,qend,qrot, zstart, zend, period, false);

if runMoments
   obj = obj.RunMoments(); 
end

end



% 0.16 (16 cm) - length of cell
% db = 0.15
% quad length & drift spacing 2x
% 7.6e-6

% .010028 - x
% .0074511 - y