function [EM,AM,NormEM,ErrEM,d] = EMMPerfect(a,d,M,w)
%DESCRIPTION: Solving electromagnetic scattering problem in 3D with M
%perfectly conducting particles
%SYNTAX     : EMMPerfect(a,d,M)
%INPUT      : a   : The radius of particles 
%             M   : Number of particles (Number of equations and unknows)
%OUTPUT     : EM     : The solution E to the scattering problem in vector form (x,y,z)
%             AM     : curl(EM)
%             NormEM : Norm of the solution EM 
%             ErrEM  : Error of the solution EM 
%             d      : Distance between neigboring bodies used 
%AUTHOR     : NHAN TRAN - nhantran@math.ksu.edu

global b alpha ES mu PI4 k k2 gamma M2 M3 M4 M5

% INITIALIZING SOME CONSTS:
% Speed of EM radio wave in free space
c = 3*10^10;
% Frequency in optics
%w = 5*10^14;
% Wave number k = 2pi/lambda
k = 2*pi*w/c;
k2 = k^2;
ik = 1i*k;
% Volume of the domain that contains all particles
VolQ = 1;
% alpha is a unit vector that indicates the direction of plane wave
alpha = [1,0,0];
% ES is E_0(0), ES \dot alpha = 0
ES = [0,1,0];

PI4 = 4*pi;
% Magnetic permebility:
mu = 1;
% Continuous distribution function of particles
N = M*(a^3)/VolQ;
% Number of particles on a side of a cube of size 1
b = round(M^(1.0/3));
% Distance between two particles: d = O(1/M^(1/3))
%d = 1/M^(1/3);
%d = 1/(b-1);
%Matrix gamma 3x3
[~,~,~,gamma,~,~,~] = EM1PerfectSphere(a,M,w,[a*10,a*10,a*10],0);
M2 = 2*M;
M3 = 3*M;
M4 = 4*M;
M5 = 5*M;
M6 = 6*M;
CD = PI4/3;

printInputs();
fprintf('\nSOLVING ELECTROMAGNETIC SCATTERING PROBLEM BY %d PERFECTLY CONDUCTING SMALL PARTICLES:\n',M);

tic
%[x,y,z] = DistributeParticles();
Pos = ParticlePosition();

%SOLVING THE PROBLEM:
A0 = SetupRHS(M6);
A = SetupMatrix(M6);
% ConditionA = cond(A)
% normA = norm(A)
X = A\A0;
%X = gmres(A,A0,10,10^-6,100);
toc
 
fprintf('\nSOLUTIONS IN (x,y,z) FOR ALL PARTICLES:\n');
sx = X(1:M,:);
sy = X(M+1:M2,:);
sz = X(M2+1:M3,:);
AM = [sx,sy,sz];
sx = X(M3+1:M4,:);
sy = X(M4+1:M5,:);
sz = X(M5+1:M6,:);
EM = [sx,sy,sz];

NormAM = norm(AM,2);
fprintf('\nNorm of curlE: %E',NormAM);

NormEM = norm(EM,2);
fprintf('\nNorm of E: %E',NormEM);

ErrEM = Error();
fprintf('\nError of E: %E\n',ErrEM);

Visualize(EM,Pos);

disp('DONE!');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function E0 = E0(PointI)
        XX = Pos(PointI,:);
        E0 = ES*exp(ik*dot(alpha,XX));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function A = SetupMatrix(MMatrixSizeM)
        A = zeros(MMatrixSizeM);                
        a32 = 2*a^3;
        CD2a3 = a32*CD;                
        
        %Run from particle 1 -> particle M:
        for j=1:M            
            for s=1:M
                if(s==j)
                    continue;
                end
                rv = Pos(j,:)-Pos(s,:);                
                r = norm(rv,2);
                r2 = r*r;
                G = exp(ik*r)/(PI4*r); %Green function
                k2G = k2*G;
                GG = ((ik*r-1)/r2)*G*rv; %GradGreen
                
                %A(j,s) = (x,y,z)
                %First row:
                row = (j-1)*6+1;
                A(row,s)=0;%For X
                A(row,s+M)=-GG(3)*CD2a3;%For Y
                A(row,s+M2)=GG(2)*CD2a3;%For Z
                A(row,s+M3)=a32*(GG(2)*gamma(3,1)-GG(3)*gamma(2,1));%For U
                A(row,s+M4)=a32*(GG(2)*gamma(3,2)-GG(3)*gamma(2,2));%For V
                A(row,s+M5)=a32*(GG(2)*gamma(3,3)-GG(3)*gamma(2,3));%For W
                
                %Second row:
                row = row+1;
                A(row,s)=GG(3)*CD2a3;%For X
                A(row,s+M)=0;%For Y
                A(row,s+M2)=-GG(1)*CD2a3;%For Z
                A(row,s+M3)=a32*(-GG(1)*gamma(3,1)-GG(3)*gamma(1,1));%For U
                A(row,s+M4)=a32*(-GG(1)*gamma(3,2)-GG(3)*gamma(1,2));%For V
                A(row,s+M5)=a32*(-GG(1)*gamma(3,3)-GG(3)*gamma(1,3));%For W
                
                %Third row:
                row = row+1;
                A(row,s)=-GG(2)*CD2a3;%For X
                A(row,s+M)=GG(1)*CD2a3;%For Y
                A(row,s+M2)=0;%For Z
                A(row,s+M3)=a32*(GG(1)*gamma(2,1)-GG(2)*gamma(1,1));%For U
                A(row,s+M4)=a32*(GG(1)*gamma(2,2)-GG(2)*gamma(1,2));%For V
                A(row,s+M5)=a32*(GG(1)*gamma(2,3)-GG(2)*gamma(1,3));%For W
                
                PGG1 = PartialGradGreen(1,r,r2,rv,G);
                PGG2 = PartialGradGreen(2,r,r2,rv,G);
                PGG3 = PartialGradGreen(3,r,r2,rv,G);
                
                %Fourth row:
                row = row+1;
                A(row,s)=CD2a3*(k2G+PGG1(1));%For X
                A(row,s+M)=CD2a3*PGG2(1);%For Y
                A(row,s+M2)=CD2a3*PGG3(1);%For Z
                A(row,s+M3)=a32*(k2G*gamma(1,1)+gamma(1,1)*PGG1(1)+gamma(2,1)*PGG2(1)+gamma(3,1)*PGG3(1));%For U
                A(row,s+M4)=a32*(k2G*gamma(1,2)+gamma(1,2)*PGG1(1)+gamma(2,2)*PGG2(1)+gamma(3,2)*PGG3(1));%For V
                A(row,s+M5)=a32*(k2G*gamma(1,3)+gamma(1,3)*PGG1(1)+gamma(2,3)*PGG2(1)+gamma(3,3)*PGG3(1));%For W 
                
                %Fifth row:
                row = row+1;
                A(row,s)=CD2a3*PGG1(2);%For X
                A(row,s+M)=CD2a3*(k2G+PGG2(2));%For Y
                A(row,s+M2)=CD2a3*PGG3(2);%For Z
                A(row,s+M3)=a32*(k2G*gamma(2,1)+gamma(1,1)*PGG1(2)+gamma(2,1)*PGG2(2)+gamma(3,1)*PGG3(2));%For U
                A(row,s+M4)=a32*(k2G*gamma(2,2)+gamma(1,2)*PGG1(2)+gamma(2,2)*PGG2(2)+gamma(3,2)*PGG3(2));%For V
                A(row,s+M5)=a32*(k2G*gamma(2,3)+gamma(1,3)*PGG1(2)+gamma(2,3)*PGG2(2)+gamma(3,3)*PGG3(2));%For W
                
                %Sixth row:
                row = row+1;
                A(row,s)=CD2a3*PGG1(3);%For X
                A(row,s+M)=CD2a3*PGG2(3);%For Y
                A(row,s+M2)=CD2a3*(k2G+PGG3(3));%For Z
                A(row,s+M3)=a32*(k2G*gamma(3,1)+gamma(1,1)*PGG1(3)+gamma(2,1)*PGG2(3)+gamma(3,1)*PGG3(3));%For U
                A(row,s+M4)=a32*(k2G*gamma(3,2)+gamma(1,2)*PGG1(3)+gamma(2,2)*PGG2(3)+gamma(3,2)*PGG3(3));%For V
                A(row,s+M5)=a32*(k2G*gamma(3,3)+gamma(1,3)*PGG1(3)+gamma(2,3)*PGG2(3)+gamma(3,3)*PGG3(3));%For W

            end       
            
            %Set the diagonal when j==s:
            row = (j-1)*6+1;
            A(row,j+M3) = 1;            
            row = row+1;
            A(row,j+M4) = 1;            
            row = row+1;
            A(row,j+M5) = 1;            
            row = row+1;
            A(row,j) = 1;            
            row = row+1;
            A(row,j+M) = 1;            
            row = row+1;
            A(row,j+M2) = 1;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function RHS = SetupRHS(MMatrixSize)
        RHS = zeros(MMatrixSize,1);     
        for s=1:M                       
            E0Vec = E0(s);
            e = ik*exp(ik*dot(alpha,Pos(s,:)));
            RHS(s) = E0Vec(1);
            RHS(s+M) = E0Vec(2);
            RHS(s+M2) = E0Vec(3);
            RHS(s+M3) = e*(ES(3)*alpha(2)-ES(2)*alpha(3));
            RHS(s+M4) = -e*(ES(3)*alpha(1)-ES(1)*alpha(3));
            RHS(s+M5) = e*(ES(2)*alpha(1)-ES(1)*alpha(2));
            break;
        end       
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Q = Qm(m)     
        Q = -2*(a^3)*(CD*AM(m,:)+(gamma*EM(m,:)')');      
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Err = Error()
        Err = (1/PI4)*(1/d^3 + k2/d + k/d^2)*a;
        sum = 0;
        for s=1:M
            sum = sum + norm(Qm(s));
        end
        
        Err = Err*sum;
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

    function G = Green(r)
        % Create a Green function in 3D       
        
        % Distance from particle s to particle t in 3D
        %r = norm(Pos(s,:)-Pos(t,:),2);
        G = exp(ik*r)/(PI4*r);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function GG = GradGreen(r,rv,G)
        %rv = Pos(PointI,:)-Pos(PointJ,:);
        %r = norm(rv,2);
        %G = exp(k1i*r)/(PI4*r);%Green function  
        
        %GG = ((k1i*G-G/r)/r)*rv;
        GG = ((ik*r-1)/r^2)*G*rv;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function PGG = PartialGradGreen(PartialIndex,r,r2,rv,G)
        % Create a Partial derivative of Grad of Green function in 3D       
        
        % Distance from particle s to particle t in 3D                
        %rv = Pos(s,:)-Pos(t,:);
        %r = norm(rv,2);   
        %G = exp(k1i*r)/(PI4*r);%Green function
        
        F = G*(ik-1/r);        
        PartialF = (-k2 - 2*ik/r + 2/(r2))*G*rv(PartialIndex)/r;        

        r0 = rv/r;
        
        Partialr0 = zeros(1,3);
        for DimIndex=1:3
            if(PartialIndex == DimIndex)
                Partialr0(DimIndex) = 1/r - (rv(PartialIndex)^2)/(r2*r);
            else
                Partialr0(DimIndex) = -(rv(DimIndex)/(r2*r))*2*rv(PartialIndex);
            end
        end
        
        PGG = PartialF*r0+F*Partialr0;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function printInputs()
        % Display the inputs arguments of wave scattering problem
        
        str = 'INPUT PARAMETERS:';
        str = strcat(str, sprintf('\nSpeed of EM wave, c: %E cm/sec',c));
        str = strcat(str, sprintf('\nFrequency, f: %E Hz',w));
        str = strcat(str, sprintf('\nWave number, k = 2pi/lambda: %E cm^-1',k));  
        str = strcat(str, sprintf('\nWave length, lambda: %E cm',2*pi/k));
        str = strcat(str, sprintf('\nVolume of the domain that contains all particles: %f cm^3',VolQ));                            
        str = strcat(str, sprintf('\nDirection of plane wave, alpha: (%s)',num2str(alpha)));
        str = strcat(str, sprintf('\nIncident field vector E_0: (%s)',num2str(ES)));
        str = strcat(str, sprintf('\nContinuous distribution function of particles: %s',num2str(N)));
        str = strcat(str, sprintf('\nMagnetic permebility: mu(x)= %s',num2str(mu)));
        str = strcat(str, sprintf('\nNumber of particles, M: %d',M));  
        str = strcat(str, sprintf('\nThe smallest distance between two neighboring particles, d: %E cm',d));
        str = strcat(str, sprintf('\nRadius of one particle, a: %E cm',a));        
        disp(str);
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Visualize(EM,Pos)
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
        %axis([0 1 0 1 0 1]);
        xlabel('x-axis');
        ylabel('y-axis');
        zlabel('z-axis');
        box on;
        axis tight;
        grid on;
        hold off;        
    end

end