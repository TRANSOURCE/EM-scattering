function [Ea,Ee,J,Ediff,M] = EM1PerfectEllipsoid(a,b,c,m,w,X,vis)
%DESCRIPTION: Solving electromagnetic wave scattering problem in 3D with
%only 1 perfectly conducting spherical body located at (0,0,0)
%SYNTAX     : [EMasym,EMexact,J,Gamma,EMdiff,Qdiff] = EM1Perfect(a,M,X,vis)
%INPUT      : a    : The radius of the particle
%             M    : Total collocation points
%             X    : A point outside the body to compute E
%             vis  : Visualize or not
%OUTPUT     : Ea     : The asymptotic electric field, i.e the solution to the EM scattering problem in vector form (x,y,z)
%             Ee     : The exact electric field
%             J      : Vector J in the solution E
%             Ediff  : Difference between the asymptotic and exact E
%             Qdiff  : Difference between the asymptotic and exact Q
%AUTHOR     : NHAN TRAN - nhantran@math.ksu.edu

global mu C k PI4 Delta

% INITIALIZING SOME CONSTS:
PI2 = 2*pi;
PI4 = 4*pi;
% Speed of EM radio wave in free space in cm
C = 3*10^10;
% Frequency in optics
if(nargin<5)
    w = 5*10^14;
end
% Wave number k = 2pi/lambda
k = PI2*w/C;
ik = 1i*k;
% characteristic constant of surface area of a ball: S=4*pi*R^2
%CS = PI4;
CD = PI4/3;
% alpha is a unit vector that indicates the direction of the incident field
alpha = [0,1,0];
% ES is E_0(0) or script_E, ES \dot alpha = 0
ES = [1,0,0];
% Constants for electric field E and magnetic field H
mu = 1;
if(nargin<6)
   X = [a,b,c]*10; 
end
MaxR = max(a,max(b,c));

fprintf('SOLVING EM SCATTERING PROBLEM BY 1 PERFECTLY CONDUCTING ELLIPSOID BODY:\n');

tic
P2TP = Point2ThetaPhi();
NsVec = NSVec();
Pos = PointPosOnBody();
Delta = ((((a*b)^1.6+(a*c)^1.6+(c*b)^1.6)/3)^(1/1.6))*PI4/(M-1);
printInputs(ES,alpha);
E0V = E0Vec(ES,alpha);
F = RHSVec(E0V);

fprintf('RESULTS:');
%fprintf('\nVector J(at M points on the spherical body) in the form (x,y,z):');
[J,err] = FindJ();
JdotN = max(checkTangential(J,NsVec));
fprintf('\nJ tangential to S(Body): \t\t (J,N) = %E',JdotN);
fprintf('\nRelative error of solving the LAS for J: %E',err);
%fprintf('\nIs J correct, relative error: \t\t\t %E',CheckErrorJ());

Gamma = [-1/3,0,0;0,-1/3,0;0,0,1/6];
Qe = Qexact(J)
% Qa = Qasym(ES,alpha)
% Qdiff = norm(Qe-Qa)/norm(Qe);
% fprintf('\nQexact vs Qasymptotic: %E',Qdiff);

RHSQ = -CD*(a*b*c)*curlE0(ES,alpha,[0,0,0]);
GQ = Xint2();
ValidateJQ = norm(Qe+GQ-RHSQ)/norm(RHSQ);
% ValidateJQ = norm(Qe+(Gamma*(Qe.')).'-RHSQ)/norm(RHSQ);
fprintf('\nValidate J and Q, relative error: \t\t %E',ValidateJQ);

%fprintf('\nVector E(at M points on the spherical body) in the form (x,y,z):');
Ea = Easym(ES,alpha,X,Qe);
Ee = Eexact(ES,alpha,X,J);
Ediff = norm(Ea-Ee)/norm(Ee);
fprintf('\nEexact(X) vs Easymp(X): \t\t\t\t %E',Ediff);
fprintf('\nDistance(X,BodyCenter): \t\t\t\t %E',norm(X));

if(nargin<7)
   vis = 0; 
end
if(vis)
    Visualize();
end

fprintf('\nDONE!\n');
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [J, error] = FindJ()
        A = MainMat();
        %S = A\F;
        %error = norm(A*S-F)/norm(F);
        [S,~,error] = gmres(A,F);        

        J = zeros(M,3);
        for i=1:M
            row = (i-1)*3+1;
            J(i,1) = S(row);
            J(i,2) = S(row+1);
            J(i,3) = S(row+2);
        end             
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function curlE0 = curlE0(ES,alpha,X)
        a0 = ik*exp(ik*dot(alpha,X));
        %curlE0 = a0*[ES(3)*alpha(2)-ES(2)*alpha(3), -ES(3)*alpha(1)+ES(1)*alpha(3), ES(2)*alpha(1)-ES(1)*alpha(2)];
        curlE0 = -a0*cross(ES,alpha);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Q = Qasym(ES,alpha)                
        BodyCenter = [0,0,0]; 
        IG1 = (eye(3)+Gamma)^-1;                
        curlE0Vec = curlE0(ES,alpha,BodyCenter).';
        Q = -CD*(a*b*c)*(IG1*curlE0Vec).';        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Q = Qexact(J)
        Q = sum(J)*Delta;
%         Q = zeros(1,3);           
%         for i=1:M   
%             Q = Q + J(i,:)*Delta;           
%         end      
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function E = Eexact(ES,alpha,X,J)
    %Compute the solution E(x) with J found from solving the linear system        
        E = zeros(1,3);
        for j=1:M
            %E = E + cross(GradGreen(X,Pos(j,:)),J(j,:))*Delta;    %cross(): VERY SLOW     
            GG = GradGreen(X,Pos(j,:));
            E = E + [GG(2)*J(j,3)-GG(3)*J(j,2),-GG(1)*J(j,3)+GG(3)*J(j,1),GG(1)*J(j,2)-GG(2)*J(j,1)];
        end
        E = E_0(ES,alpha,X)+E*Delta;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function E = Easym(ES,alpha,X,Q)
    %Compute the solution E(x) with Q found from J by solving the linear system
        BodyCenter = zeros(1,3);              
        E = E_0(ES,alpha,X) + cross(GradGreen(X,BodyCenter),Q);                   
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function E0 = E_0(ES,alpha,X)        
        E0 =  ES*exp(ik*dot(alpha,X));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function E0Vec = E0Vec(ES,alpha)
        E0Vec = zeros(M,3);        
        for i=1:M
           E0Vec(i,:) = E_0(ES,alpha,Pos(i,:)); 
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function PT = Point2ThetaPhi()   
        I = 1;        
        PT = []; 
        mp = m;        
        I1 = pi/(mp+1);
        for j=1:mp
            phi = j*I1;
            mt = floor(mp+abs(phi-pi/2)*6*mp);  
            I2 = PI2/mt;
            for i=1:mt                  
                theta = i*I2;
                PT(I,:) = [theta,phi];               
                I = I+1;
            end            
        end  
        PT =[[0,0];PT];
        PT =[PT;[0,pi]];
        M = I+1;  
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function P = Index2Point(i,j)
        P = round((i-1)*sqrt(M)+j);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [i,j] = Point2Index(I)        
        i = floor(I/sqrt(M))+1;
        j = round(I-(i-1)*sqrt(M));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The body is centered at the origin
    function P = Position(PointI)
        theta = P2TP(PointI,1);
        phi = P2TP(PointI,2);
        
        P = [a*cos(theta)*sin(phi),b*sin(theta)*sin(phi),c*cos(phi)];
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The body is centered at the origin
    function Pos = PointPosOnBody()
        Pos = zeros(M,3);        
        for i=1:M
           Pos(i,:) = Position(i); 
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function NS = NS(PointI)
        theta = P2TP(PointI,1);
        phi = P2TP(PointI,2);
        
        %NS = [b*c*cos(theta)*sin(phi),a*c*sin(theta)*sin(phi),a*b*cos(phi)]; %*sin(phi)
        NS = 2*[cos(theta)*sin(phi)/a,sin(theta)*sin(phi)/b,cos(phi)/c]; 
        NS = NS/norm(NS);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function NVec = NSVec()
        NVec = zeros(M,3);    
        for i=1:M
           NVec(i,:) = NS(i); 
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function F = RHS(E0,PointI)        
        F = 2*cross(E0(PointI,:),NsVec(PointI,:));                       
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function V = RHSVec(E0) 
        V = [];
        for i=1:M
            V = [V,RHS(E0,i)];
        end
        V = V.';
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function A = MainMat()        
        A = zeros(3*M);  
        DeltaA = 2*Delta;
        for t=1:M
            N = NsVec(t,:);  
            row = (t-1)*3+1;
            for s=1:M
                if(t==s)
                    continue;
                end
                rv = Pos(t,:)-Pos(s,:);
                r = norm(rv);
                G = exp(ik*r)/(PI4*r);    %Green function          
                GG = ((ik-1/r)*G/r)*rv;   %GradGreen                                                
                %dotp = dot(GG,N);        %VERY SLOW
                dotp = GG(1)*N(1)+GG(2)*N(2)+GG(3)*N(3);                               
                
                %A(t,s) = (Ax,Ay,Az)
                %Equation for the first coordinate Ax:                      
                col = (s-1)*3+1;
                A(row,col) = (GG(1)*N(1)-dotp)*DeltaA;  %For x
                A(row,col+1) = GG(1)*N(2)*DeltaA;       %For y
                A(row,col+2) = GG(1)*N(3)*DeltaA;       %For z
                
                %Equation for the second coordinate Ay:               
                A(row+1,col) = GG(2)*N(1)*DeltaA;
                A(row+1,col+1) = (GG(2)*N(2)-dotp)*DeltaA;
                A(row+1,col+2) = GG(2)*N(3)*DeltaA;
                
                %Equation for the third coordinate Az:                
                A(row+2,col) = GG(3)*N(1)*DeltaA;
                A(row+2,col+1) = GG(3)*N(2)*DeltaA;
                A(row+2,col+2) = (GG(3)*N(3)-dotp)*DeltaA;
            end
            %Set the diagonal when t==s:            
            A(row,row) = 1;            
            A(row+1,row+1) = 1;           
            A(row+2,row+2) = 1;                        
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function GG = GradGreen(X,Y)
        rv = X-Y;
        r = norm(rv,2);
        %G = exp(ik*r)/(PI4*r);     %Green function        
        %GG = (G*(ik-1/r)/r)*rv;
        GG = -rv/(PI4*(r*r*r));      %GG0
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function printInputs(ES,alpha)
        % Display the inputs arguments of wave scattering problem
        
        str = 'INPUT PARAMETERS:';
        str = strcat(str, sprintf('\nSpeed of EM wave, c: %0.2E cm/sec',C));
        str = strcat(str, sprintf('\nFrequency, f: %0.2E Hz',w));
        str = strcat(str, sprintf('\nWave number, k = 2pi/lambda: %0.2E cm^-1',k));             
        str = strcat(str, sprintf('\nWave length, lambda: %0.2E cm',PI2/k));                              
        str = strcat(str, sprintf('\nDirection of incident plane wave, alpha: (%s)',num2str(alpha)));
        str = strcat(str, sprintf('\nVector ES: (%s)',num2str(ES)));
        str = strcat(str, sprintf('\nMagnetic permebility: mu(x)= %s',num2str(mu)));
        str = strcat(str, sprintf('\nNumber of collocation points, M: %d',M)); 
        str = strcat(str, sprintf('\nSemi-principle axes of the body, (a,b,c): (%0.2E, %0.2E, %0.2E) cm',a,b,c));
        if(nargin > 8)                                
            %str = strcat(str, sprintf('\nBoundary impedance zeta(x)= %s',num2str(zeta)));
        end
        
        disp(str);
        fprintf('\n');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function CT = checkTangential(J,NsVec)
       CT = zeros(1,M); 
       for i=1:M
          CT(i) = norm(dot(J(i,:),NsVec(i,:))); 
       end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function LHSIE = LHSIE(s)       
        integ = zeros(1,3);         
        for t=1:M
            if(s==t)
                continue;
            end            
            rv = Pos(s,:)-Pos(t,:);
            r = norm(rv,2);            
            GG = -rv/(PI4*(r*r*r));
            integ = integ + GG*(NsVec(s,1)*J(t,1)+NsVec(s,2)*J(t,2)+NsVec(s,3)*J(t,3)) - J(t,:)*(GG(1)*NsVec(s,1)+GG(2)*NsVec(s,2)+GG(3)*NsVec(s,3));
        end               
        LHSIE = J(s,:) + 2*integ*Delta; 
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function ErrJ = CheckErrorJ()
        LHSVec = [];
        for s=1:M
            LHSVec = [LHSVec, LHSIE(s)];
        end         
        ErrJ = norm((LHSVec.')-F)/norm(F);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Xint1 = Xint1(s)
        Xint1 = zeros(1,3);   
        for t=1:M                      
            rv = Pos(s,:)-Pos(t,:);
            r = norm(rv,2);  
            if(r<=MaxR)
                continue;
            end  
            GG = -rv/(PI4*(r*r*r));             
            Xint1 = Xint1 + GG*(NsVec(s,1)*J(t,1)+NsVec(s,2)*J(t,2)+NsVec(s,3)*J(t,3));
        end               
        Xint1 = Xint1*Delta; 
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Xint2 = Xint2()
        Vec = zeros(1,3);                
        for s=1:M                              
            Vec = Vec + Xint1(s);            
        end         
        Xint2 = Vec*Delta;       
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Visualize()
        EM = zeros(M,3);
        for ii=1:M
            EM(ii,:) = Easym(ES,alpha,Pos(ii,:),Qe);
        end
        EM = real(EM);
        
        figure
        hold on;
        quiver3(Pos(:,1),Pos(:,2),Pos(:,3),EM(:,1),EM(:,2),EM(:,3));
        colormap hsv;
        view(-35,45);
        axis tight;
        xlabel('x-axis');
        ylabel('y-axis');
        zlabel('z-axis');
        grid on;
        hold off;
    end

end