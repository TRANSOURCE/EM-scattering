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

global mu c k PI4 GQVec Delta1 Delta2 GGMat F

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
GGMat = GradGreenMat();
Delta1 = pi*PI2*(a^2)/M;
Delta2 = zeros(1,M);
for ii=1:M
    [ss,tt] = Point2Index(ii);
    Delta2(ii) = Delta1*abs(sin(P2TP(ii,2)))*(ss*tt);    
end

fprintf('\nRESULT:\n');
%fprintf('\nVector J(at M points on the spherical body) in the form (x,y,z):');
[J,err] = FindJ(ES,alpha);
JdotN = max(checkTangential(J,NsVec));
fprintf('\nJ tangential to S(Body): (J,N) = %E',JdotN);
fprintf('\nIs J correct, relative error: %E',err);

Qe = Qexact(J)
GQVec = GammaQe(J);
fprintf('\nValidate J, error: %E',ValidateJ(ES,alpha,J,Qe));
Qa = Qasym(ES,alpha)
Qdiff = norm(Qe-Qa)/norm(Qe);
fprintf('\nQexact vs Qasymptotic: %E',Qdiff);

%fprintf('\nVector E(at M points on the spherical body) in the form (x,y,z):');
EMasym = Easym(ES,alpha,X,Qe);
EMexact = Eexact(ES,alpha,X,J);
EMdiff = norm(EMasym-EMexact)/norm(EMexact);
fprintf('\nEexact(X) vs Easymptotic(X): %E',EMdiff);
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

    function [J, error] = FindJ(ES,alpha)
        E0 = E0Vec(ES,alpha);
        F = RHSVec(E0);
        A = MainMat();
        S = A\F;
        error = norm(A*S-F)/norm(F);

        sx = S(1:M,:);
        sy = S(M+1:M2,:);
        sz = S(M2+1:M3,:);
        J = [sx,sy,sz];        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function curlE0 = curlE0(ES,alpha,X)
        a0 = ik*exp(ik*dot(alpha,X));
        %curlE0 = a0*[ES(3)*alpha(2)-ES(2)*alpha(3), -ES(3)*alpha(1)+ES(1)*alpha(3), ES(2)*alpha(1)-ES(1)*alpha(2)];
        curlE0 = -a0*cross(ES,alpha);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function GammaE = GammaEntry(p,q,Point)          
        GammaE = 0;
        for i=1:M                   
            GG = GGMat(i,Point,:);
            GammaE = GammaE + GG(p)*NsVec(i,q)*Delta2(i);           
        end      
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function GammaE = GammaEntryV(q,Point)          
        GammaE = zeros(1,3);
        GG = zeros(1,3);
        for i=1:M                 
            GG(1) = GGMat(i,Point,1);
            GG(2) = GGMat(i,Point,2);
            GG(3) = GGMat(i,Point,3);
            temp = NsVec(i,q)*Delta2(i);
            GammaE = GammaE + GG*temp;            
        end      
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function GQ = GammaQ(J)        
        GQ = zeros(1,3);                       
        for t=1:M
            for q=1:3
                temp = J(t,q)*Delta2(t);
                GEV = GammaEntryV(q,t);
                GQ = GQ + GEV*temp;
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function GQ = GammaQe(J)                 
        sum = zeros(1,3);
        GG = zeros(1,3);
        for s=1:M
            sum1 = zeros(1,3);
            for t=1:M               
                GG(1) = GGMat(s,t,1);
                GG(2) = GGMat(s,t,2);
                GG(3) = GGMat(s,t,3);
                dotNJ = NsVec(s,1)*J(t,1)+NsVec(s,2)*J(t,2)+NsVec(s,3)*J(t,3);
                sum1 = sum1 + GG*dotNJ*Delta2(t);
            end
            sum = sum + sum1*Delta2(s);
        end
        GQ = sum;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function RHSQ = RHSQ()
        RHSQ = zeros(1,3);
        for s=1:M
            RHSQ = RHSQ+F(s,:)*Delta2(s);
        end        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function CheckJ = ValidateJ(ES,alpha,J,Qe)       
        BodyCenter = zeros(1,3);        
        RHSQVec = -CD*(a^3)*curlE0(ES,alpha,BodyCenter);     %This might be wrong 
        normRHSQVeca = norm(RHSQVec)
        RHSQVec = RHSQ();
        normRHSQVece = norm(RHSQVec)
        CheckJ = norm(Qe+GQVec-RHSQVec)/norm(RHSQVec);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function gamma = FindGamma()        
        gamma = zeros(3);    
        sumQ = norm(Qe,2)^-2;  
        for p=1:3            
            for q=1:3
                gamma(p,q) = Qe(q)*GQVec(p)*sumQ;
            end            
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function gamma = FindGamma1()                       
        sumQ = norm(Qe,2)^-2;     
        cMin = dot(conj(Qe),GQVec)*sumQ;
        gamma = eye(3)*cMin;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function gamma = FindGamma2(Point)        
        gamma = zeros(3);        
        for p=1:3            
            for q=1:3
                gamma(p,q) = GammaEntry(p,q,Point);
            end            
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function gamma = FindGamma3()        
        gamma = -eye(3)/6;
        gamma(3,3) = -gamma(3,3);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Q = Qasym(ES,alpha)
        BodyCenter = zeros(1,3); 
        Gamma = FindGamma1();
        IG1 = (eye(3)+Gamma)^-1;
        curlE0Vec = curlE0(ES,alpha,BodyCenter)';
        Q = -CD*(a^3)*(IG1*curlE0Vec)';             
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Q = Qexact(J)
        Q = zeros(1,3);           
        for i=1:M   
            Q = Q + J(i,:)*Delta2(i);           
        end      
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function E = Eexact(ES,alpha,X,J)
    %Compute the solution E(x) with J found from solving the linear system        
        E = zeros(1,3);
        for j=1:M
            %E = E + cross(GradGreen(X,Pos(j,:)),J(j,:))*Delta2;    %cross(): VERY SLOW     
            GG = GradGreen(X,Pos(j,:));
            E = E + [GG(2)*J(j,3)-GG(3)*J(j,2),-GG(1)*J(j,3)+GG(3)*J(j,1),GG(1)*J(j,2)-GG(2)*J(j,1)]*Delta2(j);
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
        for t=1:M
            N = NsVec(t,:);            
            for col=1:M
                if(t==col)
                    continue;
                end
                
                GG = GGMat(t,col,:);
                %dotp = dot(GG,N);      VERY SLOW
                dotp = GG(1)*N(1)+GG(2)*N(2)+GG(3)*N(3);  
                DeltaA = 2*Delta2(col);                              
                
                %A(j,s) = (Ax,Ay,Az)
                %Equation for the first coordinate Ax:
                row = (t-1)*3+1;
                A(row,col) = (GG(1)*N(1)-dotp)*DeltaA;%For x
                A(row,col+M) = GG(1)*N(2)*DeltaA;%For y
                A(row,col+M2) = GG(1)*N(3)*DeltaA;%For z
                
                %Equation for the second coordinate Ay:
                row = row+1;
                A(row,col) = GG(2)*N(1)*DeltaA;
                A(row,col+M) = (GG(2)*N(2)-dotp)*DeltaA;
                A(row,col+M2) = GG(2)*N(3)*DeltaA;
                
                %Equation for the third coordinate Az:
                row = row+1;
                A(row,col) = GG(3)*N(1)*DeltaA;
                A(row,col+M) = GG(3)*N(2)*DeltaA;
                A(row,col+M2) = (GG(3)*N(3)-dotp)*DeltaA;
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

    function GGMat = GradGreenMat()
        GGMat = zeros(M,M,3);
        for s=1:M
            for t=1:M
                if(s==t)
                    continue;
                end
                GGMat(s,t,:) = GradGreen2(s,t);
            end
        end
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
        str = strcat(str, sprintf('\nVector ES: (%s)',num2str(ES)));
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