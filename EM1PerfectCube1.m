function [EMasym,EMexact,J,Gamma,EMdiff,Qdiff,Q2termSign] = EM1PerfectCube(a,M,w,X,vis)
%DESCRIPTION: Solving electromagnetic wave scattering problem in 3D with
%only 1 perfectly conducting spherical body located at (0,0,0)
%SYNTAX     : [EMasym,EMexact,J,Gamma,EMdiff,Qdiff] = EM1Perfect(a,M,X,vis)
%INPUT      : a    : The radius of the particle
%             M    : Total collocation points
%             X    : A point outside the body to compute E
%             vis  : Visualize or not
%OUTPUT     : EMasym : The asymptotic electric field, i.e the solution to the EM scattering problem in vector form (x,y,z)
%             EMexact: The exact electric field
%             J      : Vector J in the solution E
%             Gamma  : Matrix Gamma to compute Q
%             EMdiff : Difference between the asymptotic and exact E
%             Qdiff  : Difference between the asymptotic and exact Q
%AUTHOR     : NHAN TRAN - nhantran@math.ksu.edu

global mu c k cS PI4 ES alpha E0

% INITIALIZING SOME CONSTS:
PI2 = 2*pi;
PI4 = 4*pi;
% Speed of EM radio wave in free space
c = 3*10^10;
% Frequency in optics
if(nargin<3)
    w = 5*10^14;
end
% Wave number k = 2pi/lambda
k = PI2*w/c;
ik = 1i*k;
% characteristic constant of surface area of a ball: S=4*pi*R^2
cS = 4*pi;
% alpha is a unit vector that indicates the direction of the incident field
alpha = [0,1,0];
% ES is E_0(0) or script_E, ES \dot alpha = 0
ES = [1,0,0];
% Constants for electric field E and magnetic field H
mu = 1;
if(nargin<4)
   X = [a,a,a]*10; 
end

S1 = round(M/6); %First surface contains S1 collocation points
d = 2*a; %Side of the cube
m = ceil(sqrt(S1)); %Number of points on a side of the unit cube in the 1st octant, the origin is 1 of its vertex
Delta = 6*d^2/M; %Surface area of 1 small piece/subdomain
CD = 8;  %|D| = d^3 = 8*a^3

printInputs();
fprintf('SOLVING ELECTROMAGNETIC SCATTERING PROBLEM BY 1 PERFECTLY CONDUCTING PARTICLE:\n');

tic
NsVec = NSVec();
Pos = PointPosOnBody();
E0 = E0Vec();
F = RHSVec();

fprintf('\nRESULT:\n');
%fprintf('\nVector J(at M points on the spherical body) in the form (x,y,z):');
[J,err] = FindJ();
JdotN = max(checkTangential(J,NsVec));
fprintf('\nIs J tangential to the surface S of the body: (J,N) = %0.2E',JdotN);
fprintf('\nRelative error of solving the LAS for J: %E',err);

Qe = Qexact(J);
Qa = Qasym(Qe);
Qdiff = norm(Qa-Qe)/norm(Qe);
fprintf('\nQexact vs Qasymptotic: %0.2E',Qdiff);

%fprintf('\nVector E(at M points on the spherical body) in the form (x,y,z):');
EMasym = Easym(X,Qe);
EMexact = Eexact(X,J);
EMdiff = norm(EMasym-EMexact)/norm(EMexact);
fprintf('\nEexact vs Easymptotic: %0.2E',EMdiff);
fprintf('\nDistance from the point X to the center of the body: %0.2E\n',norm(X));

if(nargin<5)
   vis = 0; 
end
if(vis)
    EM3 = zeros(M,3);
    dist = 1.5;
    Y = Pos*dist;
    for ii=1:M
        EM3(ii,:) = Eexact(Y(ii,:),J);
    end
    Visualize(EM3,Y,d*dist);
end
toc

fprintf('\nDONE!\n');

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

    function curlE0 = curlE0(ES,X)
        a0 = ik*exp(ik*dot(alpha,X));
        curlE0 = a0*[ES(3)*alpha(2)-ES(2)*alpha(3), -ES(3)*alpha(1)+ES(1)*alpha(3), ES(2)*alpha(1)-ES(1)*alpha(2)];
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   function gamma = FindGamma(Q)        
        XV = [0,0,0];
        
        ESV = [1,0,0];
        RHS1 = -Q/(2*a^3)-CD*curlE0(ESV,XV);
        E0V1 = E_0(ESV,XV);
        
        ESV = [0,1,0];
        RHS2 = -Q/(2*a^3)-CD*curlE0(ESV,XV);
        E0V2 = E_0(ESV,XV);
        
        ESV = [0,0,1];
        RHS3 = -Q/(2*a^3)-CD*curlE0(ESV,XV);
        E0V3 = E_0(ESV,XV);
        
        AA = [E0V1;E0V2;E0V3];
        
        RHSV = [RHS1(1),RHS2(1),RHS3(1)].';
        R1 = AA\RHSV;
               
        RHSV = [RHS1(2),RHS2(2),RHS3(2)].';
        R2 = AA\RHSV;
        
        RHSV = [RHS1(3),RHS2(3),RHS3(3)].';
        R3 = AA\RHSV;
        
        gamma = [R1.';R2.';R3.'];        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Q = Qasym(QVec)
        BodyCenter = zeros(1,3);              
        Gamma = FindGamma(QVec)               
        
%         C = -2*(a^3); 
%         Q = C*(CD*curlE0(ES,BodyCenter)+(Gamma*E_0(ES,BodyCenter).').');         
%         Q1term = C*CD*curlE0(ES,BodyCenter);
%         Q2term = C*(Gamma*E_0(ES,BodyCenter).').';
%         Q = Q1term + Q2term;
%         Q2termSign = norm(Q2term)/(norm(Q1term)+norm(Q2term)); 
        
        IG1 = (eye(3)+Gamma)^-1;                
        curlE0Vec = curlE0(ES,BodyCenter).';
        Q = -CD*(a^3)*(IG1*curlE0Vec).'; 
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Q = Qexact(J)
        Q = zeros(1,3);           
        for i=1:M                           
            Q = Q + J(i,:)*Delta;           
        end      
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function E = Eexact(X,J)
    %Compute the solution E(x) with J found from solving the linear system        
        E = zeros(1,3);
        for j=1:M                        
            %0.2E = E + cross(GradGreen(X,Pos(j,:)),J(j,:))*Delta;    %cross(): VERY SLOW        
            GG = GradGreen(X,Pos(j,:));
            E = E + [GG(2)*J(j,3)-GG(3)*J(j,2),-GG(1)*J(j,3)+GG(3)*J(j,1),GG(1)*J(j,2)-GG(2)*J(j,1)]*Delta;
        end
        E = E_0(ES,X)+E;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function E = Easym(X,Q)
    %Compute the solution E(x) with Q found from J by solving the linear system
        BodyCenter = zeros(1,3);              
        E = E_0(ES,X) + cross(GradGreen(X,BodyCenter),Q);                   
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function E0 = E_0(ES,X)        
        E0 =  ES*exp(ik*dot(alpha,X));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function E0Vec = E0Vec()
        E0Vec = [];        
        for i=1:M
           E0Vec = [E0Vec;E_0(ES,Pos(i,:))]; 
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The body is centered at the origin
    function Pos = PointPosOnBody()
        Pos = zeros(M,3); 
        dm = d/m;
        t = 1;
        for i=1:m
            for j=1:m
                Pos(t,:) = [d,(i-0.5)*dm,(j-0.5)*dm];
                t = t+1;
            end
        end
        for i=1:m
            for j=1:m
                Pos(t,:) = [(i-0.5)*dm,d,(j-0.5)*dm];
                t = t+1;
            end
        end
        for i=1:m
            for j=1:m
                Pos(t,:) = [0,(i-0.5)*dm,(j-0.5)*dm];
                t = t+1;
            end
        end
        for i=1:m
            for j=1:m
                Pos(t,:) = [(i-0.5)*dm,0,(j-0.5)*dm];
                t = t+1;
            end
        end
        for i=1:m
            for j=1:m
                Pos(t,:) = [(i-0.5)*dm,(j-0.5)*dm,d];
                t = t+1;
            end
        end
        for i=1:m
            for j=1:m
                Pos(t,:) = [(i-0.5)*dm,(j-0.5)*dm,0];
                t = t+1;
            end
        end
        
        Pos = Pos(1:M,:); 
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function NVec = NSVec()
        NVec = zeros(M,3);          
        S2 = 2*S1;
        S3 = 3*S1;
        S4 = 4*S1;
        S5 = 5*S1;
        for i=1:S1
            NVec(i,:) = [1,0,0];
        end
        for i=S1+1:S2
            NVec(i,:) = [0,1,0];
        end
        for i=S2+1:S3
            NVec(i,:) = [-1,0,0];
        end
        for i=S3+1:S4        
            NVec(i,:) = [0,-1,0];
        end
        for i=S4+1:S5
            NVec(i,:) = [0,0,1];
        end
        for i=S5+1:M
            NVec(i,:) = [0,0,-1];
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function F = RHS(PointI)        
        F = 2*cross(E0(PointI,:),NsVec(PointI,:));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function V = RHSVec() 
        V = [];
        for i=1:M
            V = [V,RHS(i)];
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
        G = exp(ik*r)/(PI4*r);%Green function        
        GG = (G*(ik-1/r)/r)*rv;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function printInputs()
        % Display the inputs arguments of wave scattering problem
        
        str = 'INPUT PARAMETERS:';
        str = strcat(str, sprintf('\nSpeed of EM wave, c: %0.2E cm/sec',c));
        str = strcat(str, sprintf('\nFrequency, f: %0.2E Hz',w));
        str = strcat(str, sprintf('\nWave number, k = 2pi/lambda: %0.2E cm^-1',k));             
        str = strcat(str, sprintf('\nWave length, lambda: %0.2E cm',PI2/k));                              
        str = strcat(str, sprintf('\nDirection of incident plane wave, alpha: (%s)',num2str(alpha)));
        str = strcat(str, sprintf('\nVector E(0): (%s)',num2str(ES)));
        str = strcat(str, sprintf('\nMagnetic permebility: mu(x)= %s',num2str(mu)));
        str = strcat(str, sprintf('\nNumber of collocation points, M: %d',M)); 
        str = strcat(str, sprintf('\nThe side of the cube, d: %0.2E cm',d));
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

    function Visualize(EM,Pos,dis)
        EM = real(EM);
        figure
        quiver3(Pos(:,1),Pos(:,2),Pos(:,3),EM(:,1),EM(:,2),EM(:,3),2);
        hold on;        
        colormap hsv;
        view(-35,45);
        axis([0 dis 0 dis 0 dis]);
        xlabel('x-axis');
        ylabel('y-axis');
        zlabel('z-axis');
        box on;
        %axis tight;
        grid on;
        hold off;
    end

end