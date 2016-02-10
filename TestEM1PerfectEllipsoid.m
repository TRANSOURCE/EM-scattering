function TestEM1PerfectEllipsoid()
    format short;
    global w

    % Frequency in optics
    w = 5*10^14;
          
    FindE(10^-8,10^-9,10^-9,18,'EM1EllipsoidM14abc899');
    RelativeError(18,'EM1EllipsoidM14.csv');   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function FindE(a,b,c,m,SaveFile)
        fid1 = fopen(strcat(SaveFile,'_Eexact.csv'),'wt');
        fid2 = fopen(strcat(SaveFile,'_Easym.csv'),'wt');
        fid3 = fopen(strcat(SaveFile,'_Ee_Ea.csv'),'wt');
        %fid4 = fopen(strcat(SaveFile,'_J.csv'),'wt');

        fprintf(fid1,'|x-x_1|,\t E_{exact}(x),');
        fprintf(fid1,'\n');
        fprintf(fid2,'|x-x_1|,\t E_{asym}(x),');
        fprintf(fid2,'\n');
        fprintf(fid3,'|x-x_1|,\t E_{exact} vs E_{asym},');
        fprintf(fid3,'\n');
        n = 4;
        for i=1:n    
            X = [a,b,c]*10^i;
            normX = norm(X);
            [Ea,Ee,J,Ediff,M] = EM1PerfectEllipsoid(a,b,c,m,w,X,0); 

            fprintf(fid1,'%0.2E, ',normX);
            fprintf(fid1,'%s, %s, %s, ',Complex2String(Ee(1)),Complex2String(Ee(2)),Complex2String(Ee(3)));
            fprintf(fid1,'\n');

            fprintf(fid2,'%0.2E, ',normX);
            fprintf(fid2,'%s, %s, %s, ',Complex2String(Ea(1)),Complex2String(Ea(2)),Complex2String(Ea(3)));
            fprintf(fid2,'\n');

            fprintf(fid3,'%0.2E, ',normX);
            fprintf(fid3,'%0.2E, ',Ediff);
            fprintf(fid3,'\n');
        end
%         for i=1:M        
%             fprintf(fid4,'%s, %s, %s, ',Complex2String(J(i,1)),Complex2String(J(i,2)),Complex2String(J(i,3)));
%             fprintf(fid4,'\n');
%         end

        fclose(fid1);
        fclose(fid2);
        fclose(fid3);
        %fclose(fid4);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function RelativeError(m,SaveFile)
        fid = fopen(SaveFile,'wt');

        n = 4;
        a = [10^-7, 10^-8, 10^-9, 10^-10];
        b = a/10;
        c = b;
        Ediff = zeros(n,1);
        %Qdiff = zeros(n,1);        
        X = [1,1,1]*10^(-7);
        
        for i=1:n    
            [~,~,~,Ediff(i),M] = EM1PerfectEllipsoid(a(i),b(i),c(i),m,w,X,0);    
        end

        fprintf(fid,'M=%d,|x-x_1|=%0.2E, ',M,norm(X));
        fprintf(fid,'\n');

        fprintf(fid,'a, ');
        for j=1:n
            fprintf(fid,'%0.2E, ',a(j));
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'b, ');
        for j=1:n
            fprintf(fid,'%0.2E, ',b(j));
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'c, ');
        for j=1:n
            fprintf(fid,'%0.2E, ',c(j));
        end
        fprintf(fid,'\n');

        fprintf(fid,'E_{exact} vs E_{asym}, ');
        for j=1:n
            fprintf(fid,'%0.2E, ',Ediff(j));
        end
        fprintf(fid,'\n');

%         fprintf(fid,'Q_{exact} vs Q_{asym}, ');
%         for j=1:n
%             fprintf(fid,'%0.2E, ',Qdiff(j));
%         end
%         fprintf(fid,'\n');

        fclose(fid);

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