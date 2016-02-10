function [EM,NormEM,errorE,d] = EMMPerfect(a,d,M,vis)
%DESCRIPTION: Solving electromagnetic scattering problem in 3D with M
%impedance particles
%SYNTAX     : EMMImpedance(a,M)
%INPUT      : a   : The radius of particles 
%             M   : Number of particles (Number of equations and unknows)
%OUTPUT     : S   : The solution to the scattering problem in vector form (x,y,z)
%AUTHOR     : NHAN TRAN - nhantran@math.ksu.edu

global b tau w alpha ES mu PI4 k

% INITIALIZING SOME CONSTS:
% Speed of light in optics
c = 3*10^10;
% Frequency in optics
w = 5*10^14;
% Wave number k = 2pi/lambda
k = 2*pi*w/c;
ik = 1i*k;
k2 = k^2;
% Volume of the domain that contains all particles
VolQ = 1;
% alpha is a unit vector that indicates the direction of plane wave
alpha = [0,1,0];
% ES is E_0(0), ES \dot alpha = 0
ES = [1,0,0];

PI4 = 4*pi;
% Magnetic permebility:
mu = 1;
% Continuous distribution function of particles
N = M*(a^3)/VolQ;
% Number of particles on a side of a cube of size 1
b = floor(M^(1/3));
% Distance between two particles: d = O(1/M^(1/3))
%d = 1/M^(1/3);
if(nargin<3)
    d = 1/(b-1);
end

Gamma = [-1/3,0,0;0,-1/3,0;0,0,1/6];
tau = (eye(3)+Gamma)^-1;

printInputs();

fprintf('\nSOLVING ELECTROMAGNETIC SCATTERING PROBLEM BY %d SMALL PERFECTLY CONDUCTING PARTICLES:\n',M);

tic
% [x,y,z] = DistributeParticles();
Pos = ParticlePosition();

%SOLVING THE PROBLEM:
A0 = SetupRHS(M*3);
A = SetupMatrix(M*3);
% ConditionA = cond(A)
% normA = norm(A)
%X = A\A0;
[X,~,errorA] = gmres(A,A0);
fprintf('\nError of the A: %0.2E',errorA);

%fprintf('\nSOLUTIONS IN (x,y,z) FOR ALL PARTICLES:\n');
AM = zeros(M,3);
for ii=1:M
    jj = (ii-1)*3+1;
    AM(ii,:) = [X(jj),X(jj+1),X(jj+2)];
end
EM = ComputeEM(AM);
NormEM = norm(EM);
errorE = ErrorEM();
fprintf('\nNorm of the E: %0.2E',NormEM);
fprintf('\nError of the E: %0.2E\n',errorE);

if(nargin<4)
    vis = 0;
end
if(vis)
    Visualize(EM,Pos);
end

disp('DONE!');
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function E0 = E0(PointI)
        XX = Pos(PointI,:);
        %E0 = ES*exp(ik*dot(alpha,XX));        
        E0 = ES*exp(ik*(alpha(1)*XX(1)+alpha(2)*XX(2)+alpha(3)*XX(3)));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Em = ComputeEM(Am)
    %Compute the solution E(x) with Am found from solving the linear system
        Em = zeros(M,3);
        C = (PI4/3)*(a^3);
        for i=1:M
            sum = zeros(1,3);
            for j=1:M
                if(i==j)
                    continue;
                end  
                %sum = sum + cross(GradGreen(i,j),Am(j,:))*C;   %VERY SLOW
                GG = GradGreen(i,j);                
                sum = sum + [GG(2)*Am(j,3)-GG(3)*Am(j,2),-GG(1)*Am(j,3)+GG(3)*Am(j,1),GG(1)*Am(j,2)-GG(2)*Am(j,1)]*C;                
            end
            Em(i,:) = E0(i)-sum;           
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function A = SetupMatrix(MMatrixSizeM)
        A = zeros(MMatrixSizeM);
        C = (PI4/3)*(a^3); 
        C1 = C*tau(1,1);
        C2 = C*tau(2,2);
        C3 = C*tau(3,3);
        
        %Run from particle 1 -> particle M:
        for j=1:M    
            row = (j-1)*3+1;
            for s=1:M
                if(s==j)
                    continue;
                end
                rv = Pos(j,:)-Pos(s,:);
                r = norm(rv,2);
                r2 = r*r;
                G = exp(ik*r)/(PI4*r);
                PGG1 = PartialGradGreen(1,r,r2,rv,G);
                PGG2 = PartialGradGreen(2,r,r2,rv,G);
                PGG3 = PartialGradGreen(3,r,r2,rv,G);
                
                %A(j,s) = (x,y,z)
                col = (s-1)*3+1;
                %First component x:                
                A(row,col)=C1*(k2*G + PGG1(1));%For x
                A(row,col+1)=C1*PGG2(1);%For y
                A(row,col+2)=C1*PGG3(1);%For z
                
                %Second component y:                
                A(row+1,col)=C2*PGG1(2);
                A(row+1,col+1)=C2*(k2*G + PGG2(2));
                A(row+1,col+2)=C2*PGG3(2);
                
                %Third component z:                
                A(row+2,col)=C3*PGG1(3);
                A(row+2,col+1)=C3*PGG2(3);
                A(row+2,col+2)=C3*(k2*G + PGG3(3));
            end            
            %Set the diagonal when j==s:           
            A(row,row) = 1;
            A(row+1,row+1) = 1;            
            A(row+2,row+2) = 1;             
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function RHS = SetupRHS(MMatrixSize)
        RHS = zeros(MMatrixSize,1);                
        for s=1:M 
            row = (s-1)*3+1;
            e = ik*exp(ik*dot(alpha,Pos(s,:)));
            RHS(row) = tau(1,1)*e*(ES(3)*alpha(2)-ES(2)*alpha(3));
            RHS(row+1) = tau(2,2)*e*(ES(3)*alpha(1)-ES(1)*alpha(3));
            RHS(row+2) = tau(3,3)*e*(ES(2)*alpha(1)-ES(1)*alpha(2));
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Q = Qm(m)        
        Q = -(PI4/3)*(a^3)*AM(m,:);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Err = ErrorEM()
        %0.2Err = (1/PI4)*max([a/d^3 a*k2/d a*k/d^2],[],2);
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

    function GG = GradGreen(PointI,PointJ)
        rv = Pos(PointI,:)-Pos(PointJ,:);
        r = norm(rv,2);
        G = exp(ik*r)/(PI4*r);%Green function        
        GG = ((ik*G-G/r)/r)*rv;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function PGG = PartialGradGreen(PartialIndex,r,r2,rv,G)
        % Create a Partial derivative of Grad of Green function in 3D       
        
        % Distance from particle s to particle t in 3D                
        %rv = Pos(s,:)-Pos(t,:);
        %r = norm(rv,2);     
        %G = exp(ik*r)/(PI4*r);%Green function
        
        F = G*ik-G/r;        
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
        str = strcat(str, sprintf('\nSpeed of light, c: %0.2E cm/sec',c));
        str = strcat(str, sprintf('\nFrequency, f: %0.2E Hz',w));
        str = strcat(str, sprintf('\nWave number, k = 2pi/lambda: %0.2E cm^-1',k));  
        str = strcat(str, sprintf('\nWave length, lambda: %0.2E cm',2*pi/k));
        str = strcat(str, sprintf('\nVolume of the domain that contains all particles: %f cm^3',VolQ));                            
        str = strcat(str, sprintf('\nDirection of plane wave, alpha: (%s)',num2str(alpha)));
        str = strcat(str, sprintf('\nIncident field vector E_0: (%s)',num2str(ES)));
        str = strcat(str, sprintf('\nContinuous distribution function of particles: %s',num2str(N)));
        str = strcat(str, sprintf('\nNumber of particles, M: %d',M));  
        str = strcat(str, sprintf('\nThe smallest distance between two neighboring particles, d: %0.2E cm',d));
        str = strcat(str, sprintf('\nRadius of one particle, a: %0.2E cm',a));       
        disp(str);
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Visualize(EM, Pos)
        EM = real(EM);
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
        title('EM field');
        box on;
        axis tight;
        grid on;
        hold off;  
        
        figure 
        subplot(1,3,1);
        plot(EM(1:M,1));
        title('EM field: x-component');
        subplot(1,3,2);
        plot(EM(1:M,2));
        title('EM field: y-component');
        subplot(1,3,3);
        plot(EM(1:M,3));
        title('EM field: z-component');
        
        figure 
        surf(EM);
    end

end