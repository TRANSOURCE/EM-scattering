function [EMasym,EMexact,J,Gamma,EMdiff,Qdiff] = EM1Perfect(a,M,X,vis)
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

global w mu c k cS PI4 ES alpha E0

% INITIALIZING SOME CONSTS:
PI2 = 2*pi;
PI4 = 4*pi;
% Speed of EM radio wave in free space
c = 3*10^10;
% Frequency in free space
w = 10^8;
%Radio wave length in cm
%lambda = 300;
% Wave number k = 2pi/lambda
k = PI2*w/c;
ik = 1i*k;
% characteristic constant of surface area of a ball: S=4*pi*R^2
cS = 4*pi;
% alpha is a unit vector that indicates the direction of the incident field
alpha = [1,0,0];
% ES is E_0(0) or script_E, ES \dot alpha = 0
ES = [0,1,0];
% Constants for electric field E and magnetic field H
mu = 1;
rootM = sqrt(M);
M2 = 2*M;
M3 = 3*M;

printInputs(c,w,k,a,M,alpha,ES,mu);
fprintf('SOLVING ELECTROMAGNETIC SCATTERING PROBLEM BY 1 PERFECTLY CONDUCTING PARTICLE:\n');

tic
P2TP = Point2ThetaPhi();
NsVec = NSVec();
Pos = PointPosOnBody();
E0 = E0Vec();
F = RHSVec();
A = MainMat();
S = A\F;
toc

sx = S(1:M,:);
sy = S(M+1:M2,:);
sz = S(M2+1:M3,:);
fprintf('\nRESULT:\n');
%fprintf('\nVector J(at M points on the spherical body) in the form (x,y,z):');
J = [sx,sy,sz];
JdotN = max(checkTangential(J,NsVec));
fprintf('\nIs J tangential to the surface S of the body: (J,N) = %e',JdotN);

Qe = Qexact(J);
Qa = Qasym(Qe);
Qdiff = norm(Qa-Qe)/norm(Qe);
fprintf('\nThe difference between exact Q and asymptotic Q is: %e',Qdiff);

%fprintf('\nVector E(at M points on the spherical body) in the form (x,y,z):');
EMasym = Easym(X,Qe);
EMexact = Eexact(X,J);
EMdiff = norm(EMasym-EMexact)/norm(EMexact);
fprintf('\nDiference between the exact E and asymptotic E: %e',EMdiff);
fprintf('\nDistance from the point X to the center of the body: %e\n',norm(X));

tic
if(vis)
    EM3 = [];
    Y = Pos*10;
    for ii=1:M
        EM3 = [EM3;Eexact(Y(ii,:),J)];
    end
    Visualize(EM3,Y,norm(Y(1,:)));
end
toc

fprintf('\nDONE!\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function curlE0 = curlE0(X)
        a0 = ik*exp(ik*dot(alpha,X));
        curlE0 = a0*[ES(3)*alpha(2)-ES(2)*alpha(3), -ES(3)*alpha(1)+ES(1)*alpha(3), ES(2)*alpha(1)-ES(1)*alpha(2)];
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Q = Qasym(QVec)
        XX = zeros(1,3);
        C = -2*(a^3);
        QC = QVec/C;
        Gamma = [QC',QC',QC'];
        Q = C*(PI4/3*curlE0(XX)+(Gamma*ES')');           
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

    function E = Eexact(X,J)
    %Compute the solution E(x) with J found from solving the linear system
        Delta1 = pi*PI2*(a^2)/M;
        E = zeros(1,3);
        for j=1:M
            [s,t] = Point2Index(j);
            Delta2 = Delta1*abs(sin(P2TP(j,2)))*(s*t);
            E = E + cross(GradGreen(X,Pos(j,:)),J(j,:))*Delta2;    %cross(): VERY SLOW        
        end
        E = E_0(X)+E;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function E = Easym(X,Q)
    %Compute the solution E(x) with Q found from J by solving the linear system
        X1 = zeros(1,3);              
        E = E_0(X) + cross(GradGreen(X,X1),Q);                   
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function E0 = E_0(X)        
        E0 =  ES*exp(ik*dot(alpha,X));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function E0Vec = E0Vec()
        E0Vec = [];        
        for i=1:M
           E0Vec = [E0Vec;E_0(Pos(i,:))]; 
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

    function F = RHS(PointI)        
        F = 2*cross(E0(PointI,:),NsVec(PointI,:));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function V = RHSVec() 
        V = [];
        for i=1:M
            V = [V,RHS(i)];
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

    function printInputs(c,f,k,a,M,alpha,ES,mu)
        % Display the inputs arguments of wave scattering problem
        
        str = 'INPUT PARAMETERS:';
        str = strcat(str, sprintf('\nSpeed of EM wave, c: %e',c));
        str = strcat(str, sprintf('\nFrequency, f: %e',f));
        str = strcat(str, sprintf('\nWave number, k = 2pi/lambda: %e',k));             
        str = strcat(str, sprintf('\nWave length, lambda: %e',PI2/k));                              
        str = strcat(str, sprintf('\nDirection of incident plane wave, alpha: (%s)',num2str(alpha)));
        str = strcat(str, sprintf('\nVector E(0): (%s)',num2str(ES)));
        str = strcat(str, sprintf('\nMagnetic permebility: mu(x)= %s',num2str(mu)));
        str = strcat(str, sprintf('\nNumber of collocation points, M: %d',M)); 
        str = strcat(str, sprintf('\nRadius of one particle, a: %e',a));
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