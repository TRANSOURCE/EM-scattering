function S = ElecMagnetScat(a,d,M)
%DESCRIPTION: Solving electromagnetic scattering problem in 3D
%SYNTAX     : ElecMagnetScat(M)
%INPUT      : M   : Number of particles (equations and unknows)
%OUTPUT     : S   : The solution to the scattering problem in vector form (x,y,z)
%AUTHOR     : NHAN TRAN - nhantran@math.ksu.edu

% INITIALIZING SOME CONSTS:
% Speed of light in optics
c = 3*10^10;
% Frequency in optics
w = 10^14;
% Wave number k = 2pi/lambda
k = 2*pi*w/c;
% characteristic constant of surface area of a ball: S=4*pi*R^2
cS = 4*pi;
% Volume of the domain that contains all particles
Domain = 1;
% Radius of one particle
%a = 10^(-2);
% Power const with respect to the radius of particles: kappa in [0,1]
kappa = 0.9;
% Distance between two particles: d = O(a^(1/3))
Maxd = ((a^(2-kappa))/Domain)^(1/3)
% alpha is a unit vector that indicates the direction of plane wave
alpha = [1,0,0];
% ES is E_0(0), ES \dot alpha = 0
ES = [0,1,0];

PI4 = 4*pi;
% COnstants for electric field E and magnetic field H
mu = 1;
% Continuous distribution function of particles
N = M*(a^(2-kappa));
%Continuous function with Re(h) >= 0
h = 1;
% Number of particles on a side of a cube of size 1
b = ceil(M^(1/3));
% Boundary impedance
zeta = h/a^kappa;
% tau matrix
tau = 2/3;

printInputs(c,w,k,kappa,Domain,a,d,M,alpha,ES,N,h,mu,zeta);

str = 'SOLVING ELECTROMAGNETIC SCATTERING PROBLEM:';
str = strcat(sprintf('\n\n'),str);
disp(str);

tic
MatrixSize = 3*M;
[x,y,z] = DistributeParticles(d,b);
Position = ParticlePosition(M,b);

%SOLVING THE PROBLEM:
A0 = SetupRHS(MatrixSize)
A = SetupMatrix(MatrixSize)
ConditionA = cond(A)
normA = norm(A)
X = A\A0;
%X = gmres(A,A0);
toc

str = 'SOLUTIONS IN (x,y,z) FOR ALL PARTICLES:';
str = strcat(sprintf('\n\n'),str);
disp(str);
sx = X(1:M,:);
sy = X(M+1:2*M,:);
sz = X(2*M+1:MatrixSize,:);
S = [sx,sy,sz];
SolutionNorm = norm(S)
error = Error()
%disp(S);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function A = SetupMatrix(MMatrixSizeM)
        A = zeros(MMatrixSizeM);
        C = cS*(a^(2-kappa))/1i*w*mu*(tau);
        M2 = 2*M;
        %Run from particle 1 -> particle M:
        for j=1:M
            %First component x:
            row = (j-1)*3+1;
            for s=1:M
                A(row,s)=C*h*((k^2)*Green(j,s) + PartialGradGreen(1,j,s,1));
                A(row,s+M)=C*h*PartialGradGreen(2,j,s,1);
                A(row,s+M2)=C*h*PartialGradGreen(3,j,s,1);
            end
            A(row,j) = A(row,j) + 1;
            
            %Second component y:
            row = (j-1)*3+2;
            for s=1:M
                A(row,s)=C*h*PartialGradGreen(1,j,s,2);
                A(row,s+M)=C*h*((k^2)*Green(j,s) + PartialGradGreen(2,j,s,2));
                A(row,s+M2)=C*h*PartialGradGreen(3,j,s,2);
            end
            A(row,j+M) = A(row,j+M) + 1;
            
            %Third component z:
            row = (j-1)*3+3;
            for s=1:M
                A(row,s)=C*h*PartialGradGreen(1,j,s,3);
                A(row,s+M)=C*h*PartialGradGreen(2,j,s,3);
                A(row,s+M2)=C*h*((k^2)*Green(j,s) + PartialGradGreen(3,j,s,3));
            end
            A(row,j+M2) = A(row,j+M2) + 1;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function RHS = SetupRHS(MMatrixSize)
        RHS = zeros(MMatrixSize,1);
        M2 = 2*M;
        
        for s=1:M
            [xs,ys,zs] = particle2position(s,b);
            e = 1i*k*exp(1i*k*dot(alpha,[xs,ys,zs]));
            RHS(s) = e*(ES(3)*alpha(2)-ES(2)*alpha(3));
            RHS(s+M) = e*(ES(3)*alpha(1)-ES(1)*alpha(3));
            RHS(s+M2) = e*(ES(2)*alpha(1)-ES(1)*alpha(2));
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Q = Qm(m)
        Q = -zeta*((cS*a^2)/(1i*w*mu))*(tau)*S(m,:);        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Err = Error()
        Err = (1/PI4)*max([a/d^3 a*k^2/d^2],[],2);
        sum = 0;
        for s=1:M
            sum = sum + norm(Qm(s));
        end
        Err = Err*sum;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [xs,ys,zs] = particle2position(s,b)
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

    function pos = ParticlePosition(M,b)
        % Return the position in the 3D cube of particle s
        
        pos = zeros(M,3);
        for s=1:M
            [pos(s,1),pos(s,2),pos(s,3)] = particle2position(s,b);
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [x,y,z] = DistributeParticles(d,b)
        % Set the position for each particle (uniformly distributed)
        
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

    function G = Green(s,t)
        % Create a Green function in 3D
        
        if(s==t)
            G = 0;
            return;
        end
        
        % Distance from particle s to particle t in 3D
        r = sqrt((Position(s,1)-Position(t,1))^2 + (Position(s,2)-Position(t,2))^2 + (Position(s,3)-Position(t,3))^2);
        G = exp(1i*k*r)/(PI4*r);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function PGG = PartialGradGreen(PartialIndex,s,t,DimIndex)
        % Create a Partial derivative of Grad of Green function in 3D
        
        if(s==t)
            PGG = 0;
            return;
        end
        
        % Distance from particle s to particle t in 3D
        r = sqrt((Position(s,1)-Position(t,1))^2 + (Position(s,2)-Position(t,2))^2 + (Position(s,3)-Position(t,3))^2);
        
        G = Green(s,t);
        F = G*1i*k-G/r;
        PartialF = (-(k^2)*G - 2*1i*k*G/r + 2*G/(r^2))*(Position(s,PartialIndex)-Position(t,PartialIndex))/r;
        r0 = (Position(s,DimIndex)-Position(t,DimIndex))/r;
        if(PartialIndex == DimIndex)
            Partialr0 = 1/r - ((Position(s,DimIndex)-Position(t,DimIndex))^2)/(r^3);
        else
            Partialr0 = -(Position(s,DimIndex)-Position(t,DimIndex))/(r^2)*2*(Position(s,PartialIndex)-Position(t,PartialIndex));
        end
        
        PGG = PartialF*r0+F*Partialr0;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function printInputs(c,f,k,kappa,VolQ,a,d,M,alpha,ES,N,h,mu,zeta)
        % Display the inputs arguments of wave scattering problem
        
        str = 'INPUT PARAMETERS:';
        str = strcat(str, sprintf('\nSpeed of light, c: %e',c));
        str = strcat(str, sprintf('\nFrequency, f: %e',f));
        str = strcat(str, sprintf('\nWave number, k = 2pi/lambda: %e',k));
        str = strcat(str, sprintf('\nKappa: %f',kappa));
        str = strcat(str, sprintf('\nVolume of the domain that contains all particles: %f',VolQ));
        str = strcat(str, sprintf('\nRadius of one particle, a: %e',a));
        str = strcat(str, sprintf('\nDistance between two particles, d = O(a^(1/3)): %e',d));
        str = strcat(str, sprintf('\nNumber of particles, M: %d',M));      
        str = strcat(str, sprintf('\nDirection of plane wave, alpha: (%s)',num2str(alpha)));
        str = strcat(str, sprintf('\nVector E_0: (%s)',num2str(ES)));
        if(nargin > 10)
            str = strcat(str, sprintf('\nContinuous distribution function of particles: %s',num2str(N)));
            %str = strcat(str, sprintf('\nDesired refraction coefficient: %s',num2str(n)));
            %str = strcat(str, sprintf('\nOriginal refraction coefficient: %s',num2str(n0)));
            %str = strcat(str, sprintf('\nInitial field satisfies Helmholtz equation in R^3: %s',num2str(u0)));
            %str = strcat(str, sprintf('\nFunction p(x): %s',num2str(p)));
            str = strcat(str, sprintf('\nFunction h(x)= %s',num2str(h)));
            str = strcat(str, sprintf('\nFunction mu(x)= %s',num2str(mu)));
            str = strcat(str, sprintf('\nBoundary impedance zeta(x)= %s',num2str(zeta)));
        end
        
        disp(str);
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end