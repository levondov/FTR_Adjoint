% define a lattice
function obj = CreateLattice(obj, a, period, runMoments)

if nargin < 4
    runMoments = false;
end
qlength = 0.01;

qstart = [0.0, .03*a(2), .07];
db = [0.02661*a(1), -0.02661*a(3), 0.02661*a(1)];
qend = [qstart(1)+qlength, qstart(2)+qlength*2, qstart(3)+qlength];
qrot = [0.0 + 45.0*pi/180*(a(4)-1), 0.0 + 45.0*pi/180*(a(5)-1), 0.0 + 45.0*pi/180*(a(4)-1)];

obj = obj.CreateLatticeProfile(db,qstart,qend,qrot, 0, qend(3), period, false);

if runMoments
   obj = obj.RunMoments(); 
end

end