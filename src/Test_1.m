i = 1;
t = 2;
time_1 = 0;
tic
while t >= time_1 ; 
for j = 1 : 100000;

i =  i + 0.0001;

end
time_1 = toc;
end
