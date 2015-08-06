% Get last point of trajectory

coordinates = [0];

for i = 1:4
   n = nv.out(nv.wall_1.varnames{i});
   coordinates = horzcat(coordinates,n(end));
end
for i = 1:10
   n = nv.out(nv.smcec_1.varnames{i});
   coordinates = horzcat(coordinates,n(end));
end

for i = 1:4
   n = nv.out(nv.wall_2.varnames{i});
   coordinates = horzcat(coordinates,n(end));
end
for i = 1:10
   n = nv.out(nv.smcec_2.varnames{i});
   coordinates = horzcat(coordinates,n(end));
end

coordinates = coordinates(2:length(coordinates));