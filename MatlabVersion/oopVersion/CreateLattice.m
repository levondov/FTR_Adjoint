% define a lattice
function obj = CreateLattice(obj, a, period, runMoments)

if nargin < 4
    runMoments = false;
end
qlength = 0.04;

qstart = [0.0, 0.06, 0.14];
db = [0.15*a(1), -0.15*a(3), 0.15];
qend = [qstart(1)+qlength/2.0, qstart(2)+qlength, qstart(3)+qlength/2.0];
qrot = [0.0, 0.0, 0.0];

% qstart = [0.0, 0.06];
% db = [0.15*a(1), -0.15*a(3)];
% qend = [qstart(1)+qlength/2.0, qstart(2)+qlength/2.0];
% qrot = [0.0, 0.0];

obj = obj.CreateLatticeProfile(db,qstart,qend,qrot, 0, 0.16, period, false);

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