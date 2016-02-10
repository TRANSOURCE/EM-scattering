function EM = EM1smallImpedance(a,X,X1)
%DESCRIPTION: Solving electromagnetic wave scattering problem in 3D with
%only 1 impedance particle
%SYNTAX     : EM1Impedance(a,X,X1)
%INPUT      : a   : The radius of the particle
%             X   : The location in R^3 where the EM field is computed
%             X1  : Position of the particle in R^3
%OUTPUT     : E   : The electric field, i.e the solution to the EM scattering problem in vector form (x,y,z)
%AUTHOR     : NHAN TRAN - nhantran@math.ksu.edu

global zeta w mu tau c k cS PI4 ES alpha h kappa

% INITIALIZING SOME CONSTS:
% Speed of light in optics
c = 3*10^10;
% Frequency in optics
w = 10^14;
% Wave number k = 2pi/lambda
k = 2*pi*w/c;
% characteristic constant of surface area of a ball: S=4*pi*R^2
cS = 4*pi;
% Power const with respect to the radius of particles: kappa in [0,1]
kappa = 0.9;
% alpha is a unit vector that indicates the direction of plane wave
alpha = [1,0,0];
% ES is E_0(0), ES \dot alpha = 0
ES = [0,1,0];
PI4 = 4*pi;
% COnstants for electric field E and magnetic field H
mu = 1;
%Continuous function with Re(h) >= 0
h = 1;
% Boundary impedance
zeta = h/a^kappa;
% tau matrix
tau = 2/3;

VE0 = E0(X);
VQ = Q(a,X1);
VGG = GradGreen(X,X1);

% EM field:
E = VE0 + cross(VGG,VQ)

    function Q = Q(a, X)
        S = cS*a^2;
        a0 = -zeta*S/(1i*w*mu);
        Q = a0*tau*curlE0(X);
    end

    function curlE0 = curlE0(X)
        a0 = 1i*k*exp(1i*k*dot(alpha,X));
        curlE0 = a0*[ES(3)*alpha(2)-ES(2)*alpha(1), -ES(3)*alpha(1)+ES(1)*alpha(3), ES(2)*alpha(1)-ES(1)*alpha(2)];
    end

    function E0 = E0(X)
       E0 =  ES*exp(1i*k*dot(alpha,X));
    end

    function GG = GradGreen(X,X1)
        r = norm(X-X1,2);
        a0 = exp(1i*k*r);
        GG = (1i*k*a0/(PI4*r)-a0/(PI4*r^2))*(X-X1)/r;
    end
end