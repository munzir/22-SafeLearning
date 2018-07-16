%% Definitions:

% define frames:
% x0y0z0: Origin wheels-mid, z0 always up, x0 always along heading
% x1y1z1: Origin on wheels-mid, x1 along wheels-connect (L->R), y1 along
% the base (at angle -q_imu from the x0 axis)

clear all
syms x psii phi dpsi dphi dx ddpsi ddphi ddx L g mw mb real
syms XXw XYw XZw YYw YZw ZZw real
syms XXb XYb XZb YYb YZb ZZb real
syms tau_R tau_L R qimu_plus_phi qimu real
syms MXb MYb MZb real

warning off
syms t X(t) PSI(t) PHI(t) dX dPSI dPHI ddX ddPSI ddPHI real
warning on
dX=diff(X,t);dPSI=diff(PSI,t);dPHI=diff(PHI,t); 
ddX=diff(dX,t);ddPSI=diff(dPSI,t);ddPHI=diff(dPHI,t);
q = [x psii phi]'; dq = [dx dpsi dphi]'; ddq = [ddx ddpsi ddphi]';

mydiff = @(H) formula(subs(diff(symfun(subs(H,...
    [x,psii,phi,dx,dpsi,dphi,ddx,ddpsi,ddphi],...
    [X, PSI,PHI,dX,dPSI,dPHI,ddX,ddPSI,ddPHI]),t),t),...
    [X, PSI,PHI,dX,dPSI,dPHI,ddX,ddPSI,ddPHI],...
    [x,psii,phi,dx,dpsi,dphi,ddx,ddpsi,ddphi]));

q_imu = -phi + qimu_plus_phi; % qimu_plus_phi is a constant if pose fixed

%% Left Wheel


thetaL = x/R - psii*L/(2*R);
dthetaL = dx/R - dpsi*L/(2*R);
iL = [1 0 0]'; jL = [0 1 0]'; kL = [0 0 1]';
i0 = [cos(thetaL) 0 sin(thetaL)]'; j0 = [0 1 0]'; k0 = [-sin(thetaL) 0 cos(thetaL)]';
w0 = dpsi*k0;
v0 = dx*i0;
alpha0 = ddpsi*k0;
a0 = ddx*i0 + dx*(cross(dpsi*k0,i0));
rOL = (L/2)*j0;
Iw=[ZZw 0 0;0 YYw 0;0 0 ZZw];
wL = w0 + dthetaL*j0;
vGL = v0 + cross(w0, rOL);
aGL = a0 + cross(alpha0, rOL) + cross(w0, cross(w0, rOL));
HGL = Iw*wL;
p = mydiff(HGL);
dHGL = p + cross(wL,HGL);

% %% Right Wheel


thetaR = x/R + psii*L/(2*R);
dthetaR = dx/R + dpsi*L/(2*R);
iR = [1 0 0]'; jR = [0 1 0]'; kR = [0 0 1]';
i0 = [cos(thetaR) 0 sin(thetaR)]'; j0 = [0 1 0]'; k0 = [-sin(thetaR) 0 cos(thetaR)]';
w0 = dpsi*k0;
v0 = dx*i0;
alpha0 = ddpsi*k0;
a0 = ddx*i0 + dx*(cross(dpsi*k0,i0));
rOR = -(L/2)*j0;
wR = w0 + dthetaR*j0;
vGR = v0 + cross(w0, rOR);
aGR = a0 + cross(alpha0, rOR) + cross(w0, cross(w0, rOR));
HGR = Iw*wR;
dHGR = mydiff(HGR)+cross(wR,HGR);


%% Body

i1 = [1 0 0]'; j1 = [0 1 0]'; k1 = [0 0 1]';
i0 = [0 sin(q_imu) -cos(q_imu)]'; j0 = [-1 0 0]'; k0 = [0 cos(q_imu) sin(q_imu)]';
w0 = dpsi*k0;
v0 = dx*i0;
a0 = ddx*i0 + dx*(cross(dpsi*k0,i0));
IB=[XXb XYb XZb;XYb YYb YZb; XZb YZb ZZb];
wB = w0 + dphi*i1;
alphaB = ddpsi*k0 + ddphi*i1 + dphi*cross(wB, i1);
rOB = [MXb MYb MZb]'/mb;
vGB = v0 + cross(wB,rOB);
aGB = a0 + cross(alphaB, rOB) + cross(wB, cross(wB, rOB));
HGB = IB*wB;
dHGB = mydiff(HGB)+cross(wB,HGB);

%% Kanes LHS

KL = sym(zeros(3,1)); KR = sym(zeros(3,1)); KB = sym(zeros(3,1));
for i=1:3
    KL(i)=mw*aGL'*diff(vGL,dq(i))+dHGL'*diff(wL,dq(i));
    KR(i)=mw*aGR'*diff(vGR,dq(i))+dHGR'*diff(wR,dq(i));
    KB(i)=mb*aGB'*diff(vGB,dq(i))+dHGB'*diff(wB,dq(i));
end
Kw = KL + KR;
K = Kw + KB;

%% Potential Energy
gVec = [0 -g*cos(q_imu) -g*sin(q_imu)]';
V = gVec'*[MXb MYb MZb]';

%% Virtual Work

i0 = [1 0 0]'; j0 = [0 1 0]'; k0 = [0 0 1]';
dW=tau_L*j0'*(wL-dphi*j0)+tau_R*j0'*(wR-dphi*j0);

%% Frictions
syms fric_1 real

Gamma_fric = sym(zeros(3,1));
Gamma_fric(1) = -2*fric_1/R*(dq(1)/R-dq(3));
Gamma_fric(2) = -fric_1*L^2/(2*R^2)*dq(2);
Gamma_fric(3) = 2*fric_1*(dq(1)/R-dq(3));


%% Equations

AA = sym(zeros(3,3)); CC = sym(zeros(3,3)); 
QQ=sym(zeros(3,1)); Gamma=sym(zeros(3,1));
for i=1:3
    for j=1:3
        AA(i,j)=getcoeff(K(i),ddq(j),1);
        % This divides the coefficient of (dqj)(dqk) equally in all column
        % j and k 
        CC(i,j)=getcoeff(K(i), dq(j),2)*dq(j);
        ccc = getcoeff(K(i),dq(j),1); 
        CC(i,j) = CC(i,j)+ccc;
        for k=1:3
            CC(i,j) = CC(i,j) - 0.5*(getcoeff(ccc,dq(k),1))*dq(k);
        end
    end
    QQ(i) = diff(V,q(i));
    Gamma(i) = diff(dW,dq(i));
end

AA=simplify(AA);
CC=simplify(CC);

K=K-Gamma+QQ;

%% Substituting (qimu + phi) in place of qimu_plus_phi for simplification.
% We couldn't have done it before because we wanted q_imu to be
% differentiated properly because it is a function of phi

AA = subs(AA,qimu_plus_phi, qimu+phi);
CC = subs(CC,qimu_plus_phi, qimu+phi);
QQ = subs(QQ,qimu_plus_phi, qimu+phi);


%% Comparing with kim

syms mS mC Iw1 Iw2 Iw3 I1 I2 I3 d real

Acheck=simplify(subs(AA,...
    [mb,mw,L,XXw,YYw,ZZw,MXb,MYb,MZb,XXb,YYb,ZZb,XYb,XZb,YZb],...
    [mS,mC,2*L,Iw2,Iw3,Iw2,0,mS*d*cos(qimu+phi),mS*d*sin(qimu+phi),I3,...
    I2*cos(qimu)^2 + I1*sin(qimu)^2,I1*cos(qimu)^2 + I2*sin(qimu)^2,0,0,...
    I2*cos(qimu)*sin(qimu)-I1*cos(qimu)*sin(qimu)]));

Ccheck=simplify(subs(CC,...
    [mb,mw,L,XXw,YYw,ZZw,MXb,MYb,MZb,XXb,YYb,ZZb,XYb,XZb,YZb],...
    [mS,mC,2*L,Iw2,Iw3,Iw2,0,mS*d*cos(qimu+phi),mS*d*sin(qimu+phi),I3,...
    I2*cos(qimu)^2 + I1*sin(qimu)^2,I1*cos(qimu)^2 + I2*sin(qimu)^2,0,0,...
    I2*cos(qimu)*sin(qimu)-I1*cos(qimu)*sin(qimu)]));


Qcheck=simplify(subs(QQ,...
    [mb,mw,L,XXw,YYw,ZZw,MXb,MYb,MZb,XXb,YYb,ZZb,XYb,XZb,YZb],...
    [mS,mC,2*L,Iw2,Iw3,Iw2,0,mS*d*cos(qimu+phi),mS*d*sin(qimu+phi),I3,...
    I2*cos(qimu)^2 + I1*sin(qimu)^2,I1*cos(qimu)^2 + I2*sin(qimu)^2,0,0,...
    I2*cos(qimu)*sin(qimu)-I1*cos(qimu)*sin(qimu)]));

%% Final Expressions

% Acheck*[ddx ddpsi ddphi]' + Ccheck*[dx dpsi dphi]' + Qcheck = Gamma + Gamma_fric
% For definitions of symbols, refer to the paper on the link:
% https://link.springer.com/content/pdf/10.1007/s10846-005-9022-4.pdf

%% Helper Functions
function c = getcoeff(P, x, a)

[C, T] = coeffs(P, x); 
n=length(C); 
exists = 0; 
for i=1:n 
    if(isequal(T(i),x^a)); 
        exists = 1; 
        break; 
    end
end

if(exists) 
    c = C(i); 
else
    c = 0;
end
end
