function TestEMMPerfect()
Test(10,'M10.csv');
Test(125,'M125.csv');
Test(1000,'M1000.csv');
end

function Test(M,SaveFile)
fid = fopen(SaveFile,'wt');

n = 5;
a = [0.0001,0.00005,0.00001,0.000005, 0.000001];
NEM = zeros(n,1);
E = zeros(n,1);
d = zeros(n,1);
for i=1:n
    [EM,NEM(i),E(i),d(i)] = EMMImpedance1(a(i),M);   
end

fprintf(fid,'a, ');
for j=1:n
    fprintf(fid,'%0.2E, ',a(j));     
end
fprintf(fid,'\n');

fprintf(fid,'d, ');
for j=1:n
    fprintf(fid,'%0.2E, ',d(j));     
end
fprintf(fid,'\n');

fprintf(fid,'Norm of E, ');
for j=1:n
    fprintf(fid,'%0.2E, ',NEM(j));     
end
fprintf(fid,'\n');

fprintf(fid,'Error of E, ');
for j=1:n
    fprintf(fid,'%0.2E, ',E(j));     
end
fprintf(fid,'\n');

fclose(fid);

end
