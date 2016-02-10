function TestEMMPerfect()
    global w 

    % Frequency in optics
    w = 5*10^14;
    d = 10^-7;

    Test(d,27,'EMPerfectM27.csv');
    Test(d,125,'EMPerfectM125.csv');
    Test(d,1000,'EMPerfectM1000.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Test(d,M,SaveFile)
        fid = fopen(SaveFile,'wt');
        fid2 = fopen(strcat(SaveFile,'_E.csv'),'wt');

        n = 4;
        a = [10^-8,10^-9,10^-10,10^-11];
        NormEM = zeros(n,1);
        ErrEM = zeros(n,1);
        
        for i=1:n
            [EM,NormEM(i),ErrEM(i),d] = EMMPerfect(a(i),d,M);
        end
        
        fprintf(fid,'M=%d,d=%0.2E, ',M,d);
        fprintf(fid,'\n');

        fprintf(fid,'a, ');
        for j=1:n
            fprintf(fid,'%0.2E, ',a(j));
        end
        fprintf(fid,'\n');

        fprintf(fid,'a/d, ');
        for j=1:n
            fprintf(fid,'%0.2E, ',a(j)/d);
        end
        fprintf(fid,'\n');

        fprintf(fid,'Norm of E, ');
        for j=1:n
            fprintf(fid,'%0.2E, ',NormEM(j));
        end
        fprintf(fid,'\n');

        fprintf(fid,'Error of E, ');
        for j=1:n
            fprintf(fid,'%0.2E, ',ErrEM(j));
        end
        fprintf(fid,'\n');        
        
        for i=1:M        
            fprintf(fid2,'%s, %s, %s, ',Complex2String(EM(i,1)),Complex2String(EM(i,2)),Complex2String(EM(i,3)));
            fprintf(fid2,'\n');
        end
        
        fclose(fid);
        fclose(fid2);

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function s = Complex2String(C)
        re = real(C);
        im = imag(C);
        sign = '';
        if(im>=0)
            sign = '+';
        end
        s = sprintf('%0.2E%s%0.2Ei',re,sign,im); 
    end

end