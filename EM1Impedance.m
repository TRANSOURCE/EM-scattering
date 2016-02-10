function EM = EM1Impedance(a,M)
%DESCRIPTION: Solving electromagnetic wave scattering problem in 3D with
%only 1 impedance spherical body located at (0,0,0)
%SYNTAX     : EM1Impedance(a,M)
%INPUT      : a    : The radius of the particle
%             M    : Total collocation points
%OUTPUT     : E    : The electric field, i.e the solution to the EM scattering problem in vector form (x,y,z)
%AUTHOR     : NHAN TRAN - nhantran@math.ksu.edu

global zeta w mu tau c k cS PI4 ES alpha h kappa 

% INITIALIZING SOME CONSTS:
PI2 = 2*pi;
PI4 = 4*pi;
% Speed of light in optics
c = 3*10^10;
% Frequency in optics
w = 10^14;
% Wave number k = 2pi/lambda
k = PI2*w/c;
k2 = k^2;
ik = 1i*k;
% characteristic constant of surface area of a ball: S=4*pi*R^2
cS = 4*pi;
% Power const with respect to the radius of particles: kappa in [0,1]
kappa = 0.9;
% alpha is a unit vector that indicates the direction of the incident field
alpha = [1,0,0];
% ES is E_0(0) or script_E, ES \dot alpha = 0
ES = [0,1,0];
% Constants for electric field E and magnetic field H
mu = 1;
%Continuous function with Re(h) >= 0
h = 1;
% Boundary impedance
zeta = h/a^kappa;
% tau matrix
tau = 2/3;
rootM = sqrt(M);
M2 = 2*M;
M3 = 3*M;

printInputs(c,w,k,kappa,a,M,alpha,ES,mu);
fprintf('SOLVING ELECTROMAGNETIC SCATTERING PROBLEM BY 1 IMPEDANCE PARTICLE:\n');

tic
P2TP = Point2ThetaPhi();
Pos = PointPosOnBody();
NVec = NSVec();
E0 = E0Vec();
RHS = RHSVec();
AA = MainMat();
S = AA\RHS;
%S = gmres(AA,RHS);
toc

sx = S(1:M,:);
sy = S(M+1:M2,:);
sz = S(M2+1:M3,:);
%fprintf('\nRESULT:\nThe solution J(at M points on the spherical body) in the form (x,y,z):\n');
J = [sx,sy,sz]; 
JdotN = max(checkTangential(J,NVec));
fprintf('\nIs J tangential to the surface S of the body: (J,N) = %0.2E',JdotN);

Q1 = Qexact(J);
%fprintf('\nVector E(at M points on the spherical body) in the form (x,y,z):');
EM = Easym(Q1);

figure
[x, y, z] = sphere;
%f = surfl(x, y, z);
f = surf(x, y, z,'EdgeColor','none','LineStyle','none','FaceLighting','phong'); 
set(f, 'FaceAlpha', 0.3)
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

% figure
% plot3(Pos(:,1),Pos(:,2),Pos(:,3),EM(:,1),EM(:,2),EM(:,3));
% grid on;
% axis square;
% figure
% surf(real(EM));
% figure
% surf(imag(EM));

fprintf('\nDONE!\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Q = Qexact(J)
        Q = zeros(1,3);   
        Delta1 = pi*PI2*(a^2)/M;
        for i=1:M   
            [s,t] = Point2Index(i);
            %Q = Q + J(i,:)*Delta(s,t);           
            Delta2 = Delta1*abs(sin(P2TP(i,2)))*(s*t);
            Q = Q + J(i,:)*Delta2;           
        end      
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function E = Easym(Q)
    %Compute the solution E(x) with Q found from J by solving the linear system
        E = zeros(M,3);          
        for i=1:M
            E(i,:) = E0(i,:) + cross(GradGreen2(i),Q);           
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

    function [i,j] = Point2Index(I)        
        i = floor(I/rootM)+1;
        j = round(I-(i-1)*rootM);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function P = Position(PointI)
        theta = P2TP(PointI,1);
        phi = P2TP(PointI,2);
        
        P = [a*cos(theta)*sin(phi),a*sin(theta)*sin(phi),a*cos(phi)];
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Pos = PointPosOnBody()
        Pos = zeros(0);
        for i=1:M
            Pos = [Pos;Position(i)];
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function E_0 = E_0(PointI)
        X = Pos(PointI,:);
        E_0 =  ES*exp(ik*dot(alpha,X));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function E0Vec = E0Vec()
        E0Vec = [];        
        for i=1:M
           E0Vec = [E0Vec;E_0(i)]; 
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function NS = NS(PointI)
        theta = P2TP(PointI,1);
        phi = P2TP(PointI,2);        
        
        NS = [-cos(theta)*sin(phi),-sin(theta)*sin(phi),-cos(phi)];
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function NSVec = NSVec()
        NSVec = zeros(0);        
        for i=1:M
           NSVec = [NSVec;NS(i)]; 
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function curlE0 = curlE0(X)
        a0 = ik*exp(ik*dot(alpha,X));
        curlE0 = a0*[ES(3)*alpha(2)-ES(2)*alpha(3), -ES(3)*alpha(1)+ES(1)*alpha(3), ES(2)*alpha(1)-ES(1)*alpha(2)];
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function RHSV = RHSV(PointI)
        Ns = NVec(PointI,:);
        RHSV = -cross(Ns,cross(E0(PointI,:),Ns))+(zeta/1i*w*mu)*cross(Ns,curlE0(Pos(PointI,:)));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function V = RHSVec() 
        V = zeros(0);
        for i=1:M
            V = [V,RHSV(i)];
        end
        V = V';
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function D = Distance(s,t)                
        % Distance from point s to point t in 3D
        D = norm(Pos(s,:)-Pos(t,:),2);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function G = Green(r)
        % Create a Green function in 3D
        
        %r = Distance(s,t);
        G = exp(ik*r)/(PI4*r);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function AA = MainMat()
        AA = zeros(M3);
        CC = zeta/(1i*w*mu);        
        CCk2 = CC*k2;
        D = 2*(pi^2)*(a^2)/M;
        
        %Run from particle 1 -> particle M:
        for t=1:M
            N = NVec(t,:);%vectors            
            for s=1:M
                if(s==t)
                    continue;
                end
                rv = Pos(t,:)-Pos(s,:);
                r = norm(rv,2);%Distance from point t to point s 
                r2 = r^2;
                R0 = rv/r;
                G = exp(ik*r)/(PI4*r); %Green function
                CCk2G = CCk2*G;                                
                [i,j] = Point2Index(s);
                Delta = D*abs(sin(P2TP(s,2)))*(i*j);               
                                
                F = G*(ik - 1/r);
                F_r = F/r;
                FA = F*(R0(3)*N(2)-R0(2)*N(3));
                FB = F*(N(3)*R0(1)-N(1)*R0(3));
                FC = F*(R0(2)*N(1)-R0(1)*N(2));               
                
                PartialF = (G*(-k2 - 2*ik/r + 2/r2)/r)*rv;
                I1 = PartialF(1)-F*rv(1)/r2;                
                I2 = PartialF(2)-F*rv(2)/r2;
                I3 = PartialF(3)-F*rv(3)/r2;
                
                S1 = N(2)*I3;
                S2 = N(3)*I2;
                S3 = N(1)*I3;
                S4 = N(3)*I1;
                S5 = N(1)*I2;
                S6 = N(2)*I1;
                
                %AA(j,s) = (AAx,AAy,AAz)
                %0.2Equation for the first coordinate AAx:
                row = (t-1)*3+1;                
                AA(row,s)=Delta*(-N(1)*FA-CC*R0(1)*(S1-S2));%For x
                AA(row,s+M)=Delta*(-F*R0(3)-N(1)*FB+CCk2G*N(3)-CC*R0(2)*(S1-S2-N(3)*F_r));%For y
                AA(row,s+M2)=Delta*(-F*R0(2)-N(1)*FC-CCk2G*N(2)-CC*R0(3)*(S1-S2+N(2)*F_r));%For z
                
                %0.2Equation for the second coordinate AAy:
                row = row+1;
                AA(row,s)=Delta*(F*R0(3)-N(2)*FA+CCk2G*N(3)-CC*R0(1)*(-S3+S4+N(3)*F_r));
                AA(row,s+M)=Delta*(-N(2)*FB-CC*R0(2)*(-S3+S4));
                AA(row,s+M2)=Delta*(-F*R0(1)-N(2)*FC-CCk2G*N(1)-CC*R0(3)*(-S3+S4-N(1)*F_r));
                
                %0.2Equation for the third coordinate AAz:
                row = row+1;
                AA(row,s)=Delta*(-F*R0(2)-N(3)*FA+CCk2G*N(2)-CC*R0(1)*(S5-S6-N(2)*F_r));
                AA(row,s+M)=Delta*(F*R0(1)-N(3)*FB+CCk2G*N(1)-CC*R0(2)*(S5-S6+N(1)*F_r));
                AA(row,s+M2)=Delta*(-N(3)*FC-CC*R0(3)*(S5-S6));
            end
            %Set the diagonal when t==s:
            row = (t-1)*3+1;
            AA(row,t) = 1;
            row = row+1;
            AA(row,t+M) = 1;
            row = row+1;
            AA(row,t+M2) = 1; 
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function GG = GradGreen(PointI,PointJ)
        rv = Pos(PointI,:)-Pos(PointJ,:);
        r = norm(rv,2);
        G = exp(ik*r)/(PI4*r);%Green function        
        GG = ((ik*G-G/r)/r)*rv;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Return the gradient vector of Green function from Point I to the origin
    function GG = GradGreen2(PointI)
        rv = Pos(PointI,:);
        r = norm(rv,2);
        G = exp(1i*k*r)/(PI4*r);%Green function        
        GG = ((1i*k*G-G/r)/r)*rv;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function printInputs(c,f,k,kappa,a,M,alpha,ES,mu)
        % Display the inputs arguments of wave scattering problem
        
        str = 'INPUT PARAMETERS:';
        str = strcat(str, sprintf('\nSpeed of light, c: %0.2E',c));
        str = strcat(str, sprintf('\nFrequency, f: %0.2E',f));
        str = strcat(str, sprintf('\nWave number, k = 2pi/lambda: %0.2E',k));
        str = strcat(str, sprintf('\nWave length, lambda: %0.2E',PI2/k)); 
        str = strcat(str, sprintf('\nKappa: %f',kappa)); 
        str = strcat(str, sprintf('\nDirection of incident plane wave, alpha: (%s)',num2str(alpha)));
        str = strcat(str, sprintf('\nVector E(0): (%s)',num2str(ES)));        
        str = strcat(str, sprintf('\nMagnetic permebility: mu(x)= %s',num2str(mu)));
        str = strcat(str, sprintf('\nBoundary impedance zeta(x)= %s',num2str(zeta)));              
        str = strcat(str, sprintf('\nNumber of collocation points, M: %d',M));      
        str = strcat(str, sprintf('\nRadius of one particle, a: %0.2E',a));  
        
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
end