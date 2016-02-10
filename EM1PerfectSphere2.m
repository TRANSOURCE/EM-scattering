function [EMasym,EMexact,J,Gamma,EMdiff,Qdiff] = EM1PerfectSphere(a,M,w,X,vis)
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

global mu c k PI4

% INITIALIZING SOME CONSTS:
PI2 = 2*pi;
PI4 = 4*pi;
% Speed of EM radio wave in free space in cm
c = 3*10^10;
% Frequency in optics
if(nargin<3)
    w = 5*10^14;
end
% Wave number k = 2pi/lambda
k = PI2*w/c;
ik = 1i*k;
% characteristic constant of surface area of a ball: S=4*pi*R^2
%CS = PI4;
CD = PI4/3;
% alpha is a unit vector that indicates the direction of the incident field
alpha = [1,0,0];
% ES is E_0(0) or script_E, ES \dot alpha = 0
ES = [0,1,0];
% Constants for electric field E and magnetic field H
mu = 1;
if(nargin<4)
   X = [a,a,a]*10; 
end
rootM = sqrt(M);
M2 = 2*M;
M3 = 3*M;

printInputs(ES,alpha);
fprintf('SOLVING ELECTROMAGNETIC SCATTERING PROBLEM BY 1 PERFECTLY CONDUCTING PARTICLE:\n');

tic
P2TP = Point2ThetaPhi();
NsVec = NSVec();
Pos = PointPosOnBody();

fprintf('\nRESULT:\n');
%fprintf('\nVector J(at M points on the spherical body) in the form (x,y,z):');
J = FindJ(ES,alpha);
JdotN = max(checkTangential(J,NsVec));
fprintf('\nJ tangential to the surface S(Body): (J,N) = %E\n',JdotN);

Qe = Qexact(J)
Qa = Qasym(ES,alpha)
Qdiff = norm(Qa-Qe)/norm(Qe);
fprintf('\nQexact vs Qasymptotic: %E',Qdiff);
%fprintf('\nIs J correct (plug in back to the integral eq): %E',ValidateJ(ES,alpha,J,Qe));

%fprintf('\nVector E(at M points on the spherical body) in the form (x,y,z):');
EMasym = Easym(ES,alpha,X,Qe);
EMexact = Eexact(ES,alpha,X,J);
EMdiff = norm(EMasym-EMexact)/norm(EMexact);
fprintf('\nEexact vs Easymptotic: %E',EMdiff);
fprintf('\nDistance(X,BodyCenter): %E\n',norm(X));
toc

if(nargin<5)
   vis = 0; 
end
if(vis)
    EM3 = [];
    Y = Pos*1.5;
    for ii=1:M
        EM3 = [EM3;Eexact(ES,alpha,Y(ii,:),J)];
    end
    Visualize(EM3,Y,norm(Y(1,:)));
end

fprintf('\nDONE!\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function J = FindJ(ES,alpha)
        E0 = E0Vec(ES,alpha);
        F = RHSVec(E0);
        A = MainMat();
        S = A\F;

        sx = S(1:M,:);
        sy = S(M+1:M2,:);
        sz = S(M2+1:M3,:);
        J = [sx,sy,sz];        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function curlE0 = curlE0(ES,alpha,X)
        a0 = ik*exp(ik*dot(alpha,X));
        curlE0 = a0*[ES(3)*alpha(2)-ES(2)*alpha(3), -ES(3)*alpha(1)+ES(1)*alpha(3), ES(2)*alpha(1)-ES(1)*alpha(2)];
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function GammaE = GammaEntry(p,q,Point)  
        Delta1 = pi*PI2*(a^2)/M;
        GammaE = 0;
        for i=1:M   
            if(i==Point)
                continue;
            end                
            [s,t] = Point2Index(i);         
            Delta2 = Delta1*abs(sin(P2TP(i,2)))*(s*t);
            GG = GradGreen2(i,Point);
            GammaE = GammaE + GG(p)*NsVec(i,q)*Delta2;           
        end      
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function GQ = GammaQ(J)        
        GQ = zeros(1,3);       
        Delta1 = pi*PI2*(a^2)/M;           
        for t=1:M
            [m,n] = Point2Index(t);
            Delta2 = Delta1*abs(sin(P2TP(t,2)))*(m*n);
            for q=1:3
                temp = J(t,q)*Delta2;
                GQ(1) = GQ(1) + GammaEntry(1,q,t)*temp;
                GQ(2) = GQ(2) + GammaEntry(2,q,t)*temp;
                GQ(3) = GQ(3) + GammaEntry(3,q,t)*temp;
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function CheckJ = ValidateJ(ES,alpha,J,Qe)       
        BodyCenter = zeros(1,3);
        GQ = GammaQ(J);
        RHSQ = -CD*a^3*curlE0(ES,alpha,BodyCenter);
        CheckJ = norm(Qe+GQ-RHSQ)/norm(RHSQ);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function gamma = FindGamma(Point)        
        gamma = zeros(3);        
        for p=1:3            
            for q=1:3
                gamma(p,q) = GammaEntry(p,q,Point);
            end            
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function gamma = GammaMat()                
        Point = [0,0,0];
        
        alpha1 = [1,0,0];
        ES1 = [0,1,0];  
        RHS1 = -CD*(a^3)*curlE0(ES1,alpha1,Point)-Qe; %Got Qe
        
        alpha1 = [0,1,0];
        ES1 = [0,0,1];        
        J1 = FindJ(ES1,alpha1);
        Qe2 = Qexact(J1);        
        RHS2 = -CD*(a^3)*curlE0(ES1,alpha1,Point)-Qe2;
        
        alpha1 = [0,0,1];
        ES1 = [1,0,0];        
        J1 = FindJ(ES1,alpha1);
        Qe3 = Qexact(J1);        
        RHS3 = -CD*(a^3)*curlE0(ES1,alpha1,Point)-Qe3;               
        
        AA = [Qe;Qe2;Qe3];        
        %det((a^-3)*AA)
        
        RHSV = [RHS1(1),RHS2(1),RHS3(1)]';
        R1 = AA\RHSV;
               
        RHSV = [RHS1(2),RHS2(2),RHS3(2)]';
        R2 = AA\RHSV;
        
        RHSV = [RHS1(3),RHS2(3),RHS3(3)]';
        R3 = AA\RHSV;
        
        gamma = [R1';R2';R3'];
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function gamma = GammaMat0()
        %This got from Gammapq(t=0)
        C = -a*exp(ik*a)*(ik-1/a)/3;
        gamma = eye(3)*C;       
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function gamma = GammaMat1()       
        gamma = zeros(3);
        GQ = GammaQ(J);
        gamma(1,1) = GQ(1)/Qe(1);
        gamma(2,2) = GQ(2)/Qe(2);
        gamma(3,3) = GQ(3)/Qe(3);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function gamma = GammaMat2()       
        gamma = zeros(3);
        
        Point = [0,0,0];        
        alpha1 = [1,0,0];
        ES1 = [0,1,0]; 
        GQ = -CD*(a^3)*curlE0(ES1,alpha1,Point)-Qe;        
        gamma(1,1) = GQ(1)/Qe(1);
        gamma(2,2) = GQ(2)/Qe(2);
        gamma(3,3) = GQ(3)/Qe(3);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Q = Qasym0(ES,alpha)
        BodyCenter = zeros(1,3);  
        Gamma = GammaMat();       
        IG_1 = (eye(3)+Gamma)^(-1);
        Q = -CD*(a^3)*(IG_1*(curlE0(ES,alpha,BodyCenter))')';             
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Q = Qasym(ES,alpha)
        BodyCenter = zeros(1,3);  
        Gamma = GammaMat();  
        Q = -CD*(a^3)*curlE0(ES,alpha,BodyCenter)-(GammaMat*(Qe'))';             
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Q = Qexact(J)
        Q = zeros(1,3);   
        Delta1 = pi*PI2*(a^2)/M;
        for i=1:M   
            [s,t] = Point2Index(i);         
            Delta2 = Delta1*abs(sin(P2TP(i,2)))*(s*t);
            Q = Q + J(i,:)*Delta2;           
        end      
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function E = Eexact(ES,alpha,X,J)
    %Compute the solution E(x) with J found from solving the linear system
        Delta1 = pi*PI2*(a^2)/M;
        E = zeros(1,3);
        for j=1:M
            [s,t] = Point2Index(j);
            Delta2 = Delta1*abs(sin(P2TP(j,2)))*(s*t);
            %E = E + cross(GradGreen(X,Pos(j,:)),J(j,:))*Delta2;    %cross(): VERY SLOW     
            GG = GradGreen(X,Pos(j,:));
            E = E + [GG(2)*J(j,3)-GG(3)*J(j,2),-GG(1)*J(j,3)+GG(3)*J(j,1),GG(1)*J(j,2)-GG(2)*J(j,1)]*Delta2;
        end
        E = E_0(ES,alpha,X)+E;
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
        E0Vec = [];        
        for i=1:M
           E0Vec = [E0Vec;E_0(ES,alpha,Pos(i,:))]; 
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function PT = Point2ThetaPhi()
        m = ceil(rootM);
        I = 1;
        PT = zeros(M,2);
        for i=1:m
            for j=1:m
                if(I>M)
                    break;
                end                
                PT(I,1) = j*PI2/rootM; %theta
                PT(I,2) = i*pi/rootM; %phi
                I = I+1;
            end            
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function P = Index2Point(i,j)
        P = round((i-1)*rootM+j);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [i,j] = Point2Index(I)        
        i = floor(I/rootM)+1;
        j = round(I-(i-1)*rootM);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The body is centered at the origin
    function P = Position(PointI)
        theta = P2TP(PointI,1);
        phi = P2TP(PointI,2);
        
        P = a*[cos(theta)*sin(phi),sin(theta)*sin(phi),cos(phi)];
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The body is centered at the origin
    function Pos = PointPosOnBody()
        Pos = -NsVec*a;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function NS = NS(PointI)
        theta = P2TP(PointI,1);
        phi = P2TP(PointI,2);        
        
        NS = [-cos(theta)*sin(phi),-sin(theta)*sin(phi),-cos(phi)];
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function NVec = NSVec()
        NVec = [];        
        for i=1:M
           NVec = [NVec;NS(i)]; 
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
        V = V';
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function A = MainMat()        
        A = zeros(M3);         
        Delta1 = pi*PI4*(a^2)/M;
        for t=1:M
            N = NsVec(t,:);            
            for col=1:M
                if(t==col)
                    continue;
                end
                
                GG = GradGreen(Pos(t,:),Pos(col,:));
                %dotp = dot(GG,N);      VERY SLOW
                dotp = GG(1)*N(1)+GG(2)*N(2)+GG(3)*N(3);
                [i,j] = Point2Index(col);
                Delta2 = Delta1*abs(sin(P2TP(col,2)))*(i*j);                              
                
                %A(j,s) = (Ax,Ay,Az)
                %Equation for the first coordinate Ax:
                row = (t-1)*3+1;
                A(row,col) = (GG(1)*N(1)-dotp)*Delta2;%For x
                A(row,col+M) = GG(1)*N(2)*Delta2;%For y
                A(row,col+M2) = GG(1)*N(3)*Delta2;%For z
                
                %Equation for the second coordinate Ay:
                row = row+1;
                A(row,col) = GG(2)*N(1)*Delta2;
                A(row,col+M) = (GG(2)*N(2)-dotp)*Delta2;
                A(row,col+M2) = GG(2)*N(3)*Delta2;
                
                %Equation for the third coordinate Az:
                row = row+1;
                A(row,col) = GG(3)*N(1)*Delta2;
                A(row,col+M) = GG(3)*N(2)*Delta2;
                A(row,col+M2) = (GG(3)*N(3)-dotp)*Delta2;
            end
            %Set the diagonal when t==col:
            row = (t-1)*3+1;
            A(row,t) = 1;
            row = row+1;
            A(row,t+M) = 1;
            row = row+1;
            A(row,t+M2) = 1;                        
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

    function GG = GradGreen2(PointI,PointJ)
        rv = Pos(PointI,:)-Pos(PointJ,:);
        r = norm(rv,2);
        G = exp(ik*r)/(PI4*r);%Green function  
        
        %GG = ((k1i*G-G/r)/r)*rv;
        GG = ((ik*r-1)/r^2)*G*rv;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function printInputs(ES,alpha)
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
        str = strcat(str, sprintf('\nRadius of one particle, a: %0.2E cm',a));
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

    function Visualize(EM,Pos,a)
        figure
        [x, y, z] = sphere;
        % f = surfl(x, y, z);
        f = surf(x, y, z,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
        set(f, 'FaceAlpha', 0.3);
        hold on;
        quiver3(Pos(:,1),Pos(:,2),Pos(:,3),EM(:,1),EM(:,2),EM(:,3));
        colormap hsv;
        view(-35,45);
        axis([-a a -a a -a a]);
        xlabel('x-axis');
        ylabel('y-axis');
        zlabel('z-axis');
        grid on;
        hold off;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute the surface area of S_ij: 1 <= i,j <= rootM
    function D = Delta(i,j)
       D = ((a^2)*PI2*i/rootM)*abs(cos(j*pi/rootM)-cos((j+1)*pi/rootM));
    end

end