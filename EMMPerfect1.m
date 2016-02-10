function [EM,S,NEM,E,d] = EMMPerfect(a,M)
%DESCRIPTION: Solving electromagnetic scattering problem in 3D with M
%perfectly conducting particles
%SYNTAX     : EMMPerfect(a,d,M)
%INPUT      : a   : The radius of particles 
%             M   : Number of particles (Number of equations and unknows)
%OUTPUT     : EM  : The solution E to the scattering problem in vector form (x,y,z)
%             NEM : Norm of the solution EM 
%             E   : Error of the solution EM 
%             S   : Vector curl(E) 
%AUTHOR     : NHAN TRAN - nhantran@math.ksu.edu

global b w alpha ES mu PI4 k

% INITIALIZING SOME CONSTS:
% Speed of light in optics
c = 3*10^10;
% Frequency in optics
w = 10^14;
% Wave number k = 2pi/lambda
k = 2*pi*w/c;
% Volume of the domain that contains all particles
Domain = 1;
% alpha is a unit vector that indicates the direction of plane wave
alpha = [1,0,0];
% ES is E_0(0), ES \dot alpha = 0
ES = [0,1,0];

PI4 = 4*pi;
% Magnetic permebility:
mu = 1;
% Continuous distribution function of particles
N = M*(a^3)/Domain;
% Number of particles on a side of a cube of size 1
b = round(M^(1.0/3));
% Distance between two particles: d = O(1/M^(1/3))
%d = 1/M^(1/3);
d = 1/(b-1);

printInputs(c,w,k,Domain,a,d,M,alpha,ES,N,mu);
fprintf('\nSOLVING ELECTROMAGNETIC SCATTERING PROBLEM BY %d PERFECTLY CONDUCTING SMALL PARTICLES:\n',M);

tic
MatrixSize = 3*M;
%[x,y,z] = DistributeParticles();
Pos = ParticlePosition();

%SOLVING THE PROBLEM:
A0 = SetupRHS(MatrixSize);
A = SetupMatrix(MatrixSize);
% ConditionA = cond(A)
% normA = norm(A)
X = A\A0;
%X = gmres(A,A0);
toc

fprintf('\nSOLUTIONS IN (x,y,z) FOR ALL PARTICLES:\n');
sx = X(1:M,:);
sy = X(M+1:2*M,:);
sz = X(2*M+1:MatrixSize,:);
S = [sx,sy,sz];
EM = Em(S);
NEM = norm(EM);
E = Error();
fprintf('\nNorm of the solution: %e',NEM);
fprintf('\nError of the solution: %e\n',E);
%disp(S);

figure
% [x, y, z] = sphere;
% h = surfl(x, y, z); 
% set(h, 'FaceAlpha', 0.1)
quiver3(Pos(:,1),Pos(:,2),Pos(:,3),EM(:,1),EM(:,2),EM(:,3)); 
hold on;
colormap hsv;
view(-35,45);
% sd = (Domain^(1/3))/2;
% axis([-sd sd -sd sd -sd sd]);
axis([0 1 0 1 0 1]);
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis');
box on;
axis tight;
grid on;
hold off;

disp('DONE!');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function E0 = E0(PointI)
        XX = Pos(PointI,:);
        E0 =  ES*exp(1i*k*dot(alpha,XX));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Em = Em(Am)
    %Compute the solution E(x) with Am found from solving the linear system
        Em = zeros(M,3);
        C = (PI4/3)*(a^3);
        for i=1:M
            sum = zeros(1,3);
            for j=1:M
                if(i==j)
                    continue;
                end                
                sum = sum + cross(GradGreen(i,j),Am(j,:))*C;
            end
            Em(i,:) = E0(i)-sum;           
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function A = SetupMatrix(MMatrixSizeM)
        A = zeros(MMatrixSizeM);
        C = (PI4/3)*(a^3);
        M2 = 2*M;
        k2 = k^2;
        
        %Run from particle 1 -> particle M:
        for j=1:M            
            for s=1:M
                if(s==j)
                    continue;
                end
                r = norm(Pos(s,:)-Pos(j,:),2);
                G = Green(j,s,r);
                
                %A(j,s) = (x,y,z)
                %First component x:
                row = (j-1)*3+1;
                A(row,s)=C*(k2*G + PartialGradGreen(1,j,s,1,r,G));%For x
                A(row,s+M)=C*PartialGradGreen(2,j,s,1,r,G);%For y
                A(row,s+M2)=C*PartialGradGreen(3,j,s,1,r,G);%For z
                
                %Second component y:
                row = (j-1)*3+2;
                A(row,s)=C*PartialGradGreen(1,j,s,2,r,G);
                A(row,s+M)=C*(k2*G + PartialGradGreen(2,j,s,2,r,G));
                A(row,s+M2)=C*PartialGradGreen(3,j,s,2,r,G);
                
                %Third component z:
                row = (j-1)*3+3;
                A(row,s)=C*PartialGradGreen(1,j,s,3,r,G);
                A(row,s+M)=C*PartialGradGreen(2,j,s,3,r,G);
                A(row,s+M2)=C*(k2*G + PartialGradGreen(3,j,s,3,r,G));
            end            
            %Set the diagonal when j==s:
            row = (j-1)*3+1;
            A(row,j) = 1;
            row = (j-1)*3+2;
            A(row,j+M) = 1;
            row = (j-1)*3+3;
            A(row,j+M2) = 1;             
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function RHS = SetupRHS(MMatrixSize)
        RHS = zeros(MMatrixSize,1);
        M2 = 2*M;
        
        for s=1:M                       
            e = 1i*k*exp(1i*k*dot(alpha,Pos(s,:)));
            RHS(s) = e*(ES(3)*alpha(2)-ES(2)*alpha(3));
            RHS(s+M) = e*(ES(3)*alpha(1)-ES(1)*alpha(3));
            RHS(s+M2) = e*(ES(2)*alpha(1)-ES(1)*alpha(2));
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Q = Qm(m)
        D = (PI4/3)*(a^3);
        Q = -S(m,:)*D;       
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Err = Error()
        Err = (1/PI4)*(1/d^3 + k^2/d + k/d^2)*a;
        sum = 0;
        for s=1:M
            sum = sum + norm(Qm(s));
        end
        
        Err = Err*sum;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [xs,ys,zs] = particle2position(s)
        % Return the position in the 3D cube of particle s
        
        % The first particle [x1,y1,z1] is at the left end bottom corner of the
        % cube and is called particle number 1. The next one will be on the same
        % row, go to the right. When finishing the first line, go to the second line
        % and start at the first column again. When finishing the first plane, move
        % up.
        
        % % Number of particles on a side of a cube of size 1
        % b = ceil(M^(1/3));
        
        % Find the plane where the particle s is on
        plane = floor((s-1)/(b^2));
        % [x1,x2,x3] is an array index
        x3 = plane + 1;
        x2 = mod((s-1), b) + 1;
        t = mod(s-1,b^2);
        x1 = floor(t/b) + 1;
        % Find the position of [x1,x2,x3] in Cartesian coordinates
        xs = x(x1);
        ys = y(x2);
        zs = z(x3);
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [x,y,z] = DistributeParticles()
        % Set the position for each particle (uniformly distributed) in a
        % cube centered at the origin
        
        % % Number of particles on a side of a cube of size 1
        % b = ceil(M^(1/3));
        
        x0 = -0.5-d;
        y0 = -0.5-d;
        z0 = -0.5-d;
        
        % The first particle [x1,y1,z1] is at the left end bottom corner of the
        % cube and is called particle number 1.
        x = zeros(1,b);
        y = zeros(1,b);
        z = zeros(1,b);
        
        for s=1:b
            x(s) = x0 + d*s;
            y(s) = y0 + d*s;
            z(s) = z0 + d*s;
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function m = Index2Order(m1,m2,m3)
        m = (m1-1)+(m2-1)*b+(m3-1)*b^2+1;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [m1,m2,m3] = Order2Index(m)
        m3 = floor((m-1)/b^2)+1;
        red1 = mod(m-1,b^2);
        m2 = floor(red1/b)+1;
        m1 = mod(red1,b)+1;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [x,y,z] = Particle2Position(m)
        % Set the position for each particle (uniformly distributed)
        % The first particle [x1,y1,z1] is at the left end bottom corner of the
        % cube and is called particle number 1.
        % The cube is in the first octant and the origin is one of its
        % vertex
        [m1,m2,m3] = Order2Index(m);
        x = (m1-1)*d;
        y = (m2-1)*d;
        z = (m3-1)*d;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function pos = ParticlePosition()
        % Return the position in the 3D cube of particles
        
        pos = zeros(M,3);
        for s=1:M
            [pos(s,1),pos(s,2),pos(s,3)] = Particle2Position(s);
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function G = Green(s,t,r)
        % Create a Green function in 3D       
        
        % Distance from particle s to particle t in 3D
        %r = norm(Pos(s,:)-Pos(t,:),2);
        G = exp(1i*k*r)/(PI4*r);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function GG = GradGreen(PointI,PointJ)
        rv = Pos(PointI,:)-Pos(PointJ,:);
        r = norm(rv,2);
        G = exp(1i*k*r)/(PI4*r);%Green function        
        GG = ((1i*k*G-G/r)/r)*rv;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function PGG = PartialGradGreen(PartialIndex,s,t,DimIndex,r,G)
        % Create a Partial derivative of Grad of Green function in 3D       
        
        % Distance from particle s to particle t in 3D
        %r = norm(Pos(s,:)-Pos(t,:),2);
        
        %G = exp(1i*k*r)/(PI4*r);%Green function
        F = G*1i*k-G/r;
        V0 = Pos(s,PartialIndex)-Pos(t,PartialIndex);
        PartialF = (-(k^2) - 2*1i*k/r + 2/(r^2))*G*V0/r;
        V = Pos(s,DimIndex)-Pos(t,DimIndex);
        r0 = V/r;
        if(PartialIndex == DimIndex)
            Partialr0 = 1/r - (V^2)/(r^3);
        else
            Partialr0 = -(V/(r^3))*2*V0;
        end
        
        PGG = PartialF*r0+F*Partialr0;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function printInputs(c,f,k,VolQ,a,d,M,alpha,ES,N,mu)
        % Display the inputs arguments of wave scattering problem
        
        str = 'INPUT PARAMETERS:';
        str = strcat(str, sprintf('\nSpeed of light, c: %e',c));
        str = strcat(str, sprintf('\nFrequency, f: %e',f));
        str = strcat(str, sprintf('\nWave number, k = 2pi/lambda: %e',k));      
        str = strcat(str, sprintf('\nVolume of the domain that contains all particles: %f',VolQ));
        str = strcat(str, sprintf('\nRadius of one particle, a: %e',a));
        str = strcat(str, sprintf('\nThe smallest distance between two neighboring particles, d = O(1/M^(1/3)): %e',d));
        str = strcat(str, sprintf('\nNumber of particles, M: %d',M));      
        str = strcat(str, sprintf('\nDirection of plane wave, alpha: (%s)',num2str(alpha)));
        str = strcat(str, sprintf('\nIncident field vector E_0: (%s)',num2str(ES)));
        if(nargin > 10)
            str = strcat(str, sprintf('\nContinuous distribution function of particles: %s',num2str(N)));
            %str = strcat(str, sprintf('\nDesired refraction coefficient: %s',num2str(n)));
            %str = strcat(str, sprintf('\nOriginal refraction coefficient: %s',num2str(n0)));
            %str = strcat(str, sprintf('\nInitial field satisfies Helmholtz equation in R^3: %s',num2str(u0)));
            %str = strcat(str, sprintf('\nFunction p(x): %s',num2str(p)));          
            str = strcat(str, sprintf('\nMagnetic permebility: mu(x)= %s',num2str(mu)));
            %str = strcat(str, sprintf('\nBoundary impedance zeta(x)= %s',num2str(zeta)));
        end
        
        disp(str);
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end