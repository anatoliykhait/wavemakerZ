% Second-order wave generation
% Anatoliy Khait (haitanatoliy@gmail.com)
% Lev Shemer (shemerl@tauex.tau.ac.il ; shemer@eng.tau.ac.il)
%
% If you use this script please cite our papers:
% [1] A. Khait, L. Shemer.
%     Nonlinear wave generation by a wavemaker in deep to intermediate
%     water depth.
%     Ocean Engineering 182 (2019) 222–234
% [2] Anatoliy Khait, Lev Shemer.
%     Nonlinear generation of narrow-banded wave trains.
%     OMAE2019-95364

clear all
clearvars
clearvars -global
format long

global grav;
grav = 9.81;
zerovalue = 1e-100;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the wavemaker shape
%
%                       /  free surface
% ---------|-------------------------------
%     ^               /
%     |    |         /
%   h |             /
%     |    |       /
%     v           /          bottom
% ---------------/======================
%     ^    |    /
%     |        /
%  lh |    |  /
%     |      /
%     v    |/
% ---------o    Location of the hinge
%{
h = 0.6;
lh = 1e5;
%}
h = 0.74;
lh = 1e5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Definition of the target wave train at x=0
%{
T0 = 2.8;
ka = 0.06;
w0 = 2*pi/T0;
Omega = w0/10;  % w0 +- Omega

syms kk;
disp_rel = w0^2 == grav * kk * tanh(kk * h);
k0 = vpasolve(disp_rel, kk, [0 inf]);
k0 = double(k0);

Tenv = 2*pi/Omega;
a0 = ka / k0;

j = 1;
ti = 0;
dt = Tenv/512;
while ti < Tenv
    ta(j) = ti;
    Ea(j) = a0 * cos(Omega*(ti-Tenv/4)) * cos(w0*(ti-Tenv/4));
    
    ti = ti + dt;
    j  = j + 1;
end
%}
T0 = 3.0;
ka = 0.06;
w0 = 2*pi/T0;

syms kk;
disp_rel = w0^2 == grav * kk * tanh(kk * h);
k0 = vpasolve(disp_rel, kk, [0 inf]);
k0 = double(k0);

a0 = ka / k0;

j = 1;
ti = 0;
dt = T0/128;
while ti < T0*10
    ta(j) = ti;
    Ea(j) = a0 * cos(w0*ti-pi/2);
    
    ti = ti + dt;
    j  = j + 1;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ETA = fft(Ea);
ETA = ETA(1 : fix(length(ETA) / 2));
ETA = ETA / length(ETA);

speccut = 30;
ETA = ETA(2:speccut+1);
w = [1 : speccut] * 2 * pi / (ta(end)-ta(1));

syms kk;
disp('Start evaluation of k and k2w');
k = zeros(size(w));
for j = 1:length(w)
    disp_rel = w(j)^2 == grav * kk * tanh(kk * h);
    wave_num = vpasolve(disp_rel, kk, [0 inf]);
    k(j) = double(wave_num);
end
k2w = zeros(size(w));
for j = 1:length(w)
    disp_rel = (2*w(j))^2 == grav * kk * tanh(kk * h);
    wave_num = vpasolve(disp_rel, kk, [0 inf]);
    k2w(j) = double(wave_num);
end
disp('Complete evaluation of k and k2w');

% Avoid all unnecessary frequencies
for j = 1:length(ETA)
    if abs(ETA(j)) < 1e-3
        ETA(j) = zerovalue;
    end
end

% Change sign of w
% From fft we have: c exp(i w t)
% But we need: c exp(i (kx - wt)) = c (i kx - i wt)
% Therefore c should be complex conjugated:
% (a + ib) exp(i wt) = (a - ib) exp(-i wt)
Aa = conj(ETA);
clear ETA;

Ea = Ea*0;
for j = 1:length(w)
    Ea = Ea + real(Aa(j) * exp(-1i*w(j)*ta));
end

%%

% Wavemaker velocities at 1st order
kappa1 = grav*k./w./cosh(k*h);
Lambda1 = kappa1.*k*(h+lh).*(sinh(2*k*h)+2*k*h) ./ ...
          (4 * (1 - cosh(k*h) + k*(h+lh).*sinh(k*h)));
AdX1 = Lambda1 .* Aa;

dX1 = zeros(size(Ea));
for l = 1:length(w)
    dX1 = dX1 + real(AdX1(l) * exp(-1i*w(l)*ta));
end

% Wavemaker displacements at the 1st order
spec_X1 = 1i * AdX1 ./ w;

X1 = zeros(size(dX1));
for l = 1:length(w)
    X1 = X1 + real(spec_X1(l) * exp(-1i*w(l)*ta));
end

% 2nd order correction due to finite displacements of the wavemaker
AdX2 = zeros(size(AdX1));
AdX2 = AdX2 + zerovalue;
for l = 1:length(w)
    if 2*l <= speccut
        kappa_d = Lambda1(l) * grav*k(l) / (2*(h+lh)*w(l)^2 * cosh(k(l)*h));
        G = k2w(l)/(k(l)-k2w(l))^2 - k2w(l)/(k(l)+k2w(l))^2 + ...
            k2w(l)*cosh((k2w(l)+k(l))*h)/(k(l)+k2w(l))^2 + ...
            k(l)*(h+lh)*sinh((k(l)+k2w(l))*h)/(k(l)+k2w(l)) + ...
            k(l)*(k(l)-k2w(l))*(h+lh)*sinh((k(l)-k2w(l))*h)/(k(l)-k2w(l))^2 - ...
            k2w(l)*cosh((k(l)-k2w(l))*h)/(k(l)-k2w(l))^2;
        AdX2(2*l) = -kappa_d*k2w(l)^2*(h+lh) / ...
                    (2*(1-cosh(k2w(l)*h)+k2w(l)*(h+lh)*sinh(k2w(l)*h))) * G * Aa(l)^2;
    end
end

dX2 = zeros(size(dX1));
for l = 1:length(w)
    dX2 = dX2 + real(AdX2(l) * exp(-1i*w(l)*ta));
end

spec_X2 = 1i * AdX2 ./ w;

X2 = zeros(size(dX2));
for l = 1:length(w)
    X2 = X2 + real(spec_X2(l) * exp(-1i*w(l)*ta));
end

%%

% Correction due to bound waves
coord = 0;
Eb = zeros(size(Ea));
l = 1;
for j = 1:length(w)
    for m = 1:length(w)
        ki = k(j) + k(m);
        wi = sqrt(grav * ki * tanh(ki * h));
        bnd = -V11(wi,w(j),w(m),ki,k(j),k(m)) / (wi-w(j)-w(m)) * ...
               (pi * sqrt(2*grav/w(j)) * Aa(j)) * (pi * sqrt(2*grav/w(m)) * Aa(m));
        AbZ(l) = (1/pi) * sqrt(wi/(2*grav)) * bnd;
        wb(l) = w(j) + w(m);
        kb(l) = ki;
        hnum(l) = j + m;
        Eb = Eb + real(AbZ(l) * exp(1i*(kb(l)*coord - wb(l)*ta)));
        l = l + 1;
        
        if (j~=m)
            ki = - k(j) + k(m);
            wi = sqrt(grav * ki * tanh(ki * h));
            bnd = -V22(wi,w(j),w(m),ki,k(j),k(m)) / (wi+w(j)-w(m)) * ...
                   conj(pi * sqrt(2*grav/w(j)) * Aa(j)) * (pi * sqrt(2*grav/w(m)) * Aa(m));
            AbZ(l) = (1/pi) * sqrt(wi/(2*grav)) * bnd;
            wb(l) = - w(j) + w(m);
            kb(l) = ki;
            hnum(l) = - j + m;
            Eb = Eb + real(AbZ(l) * exp(1i*(kb(l)*coord - wb(l)*ta)));
            l = l + 1;
        end
        
        ki = - k(j) - k(m);
        wi = sqrt(grav * ki * tanh(ki * h));
        bnd = -V33(wi,w(j),w(m),ki,k(j),k(m)) / (wi+w(j)+w(m)) * ...
               conj(pi * sqrt(2*grav/w(j)) * Aa(j)) * conj(pi * sqrt(2*grav/w(m)) * Aa(m));
        AbZ(l) = (1/pi) * sqrt(wi/(2*grav)) * bnd;
        wb(l) = - w(j) - w(m);
        kb(l) = ki;
        hnum(l) = - j - m;
        Eb = Eb + real(AbZ(l) * exp(1i*(kb(l)*coord - wb(l)*ta)));
        l = l + 1;
    end
end

% avemaker motion to produce bound waves
AdX2bZ = zeros(size(AbZ));
spec_X2bZ = zeros(size(AbZ));
for l = 1:length(wb)
    kappa2 = 3*wb(l) / (sinh(kb(l)*h) * (2 + cosh(kb(l)*h)));
    Lambda2 = kappa2*kb(l)*(h+lh)*(sinh(2*kb(l)*h)+2*kb(l)*h) / ...
             (4 * (1 - cosh(kb(l)*h) + kb(l)*(h+lh)*sinh(kb(l)*h)));
    if ~isnan(Lambda2)
        AdX2bZ(l) = Lambda2 * AbZ(l);
        spec_X2bZ(l) = 1i * AdX2bZ(l) / wb(l);
    end
end

% Convert spectra to positive-frequecy representation
Ab = zeros(size(Aa));
AdX2b = zeros(size(AdX2));
spec_X2b = zeros(size(spec_X2));
for j = 1:length(hnum)
    if abs(hnum(j)) <= speccut
        if hnum(j) >= 0
            Ab(abs(hnum(j))) = Ab(abs(hnum(j))) + AbZ(j);
            AdX2b(abs(hnum(j))) = AdX2b(abs(hnum(j))) + AdX2bZ(j);
            spec_X2b(abs(hnum(j))) = spec_X2b(abs(hnum(j))) + spec_X2bZ(j);
        else
            Ab(abs(hnum(j))) = Ab(abs(hnum(j))) + conj(AbZ(j));
            AdX2b(abs(hnum(j))) = AdX2b(abs(hnum(j))) + conj(AdX2bZ(j));
            spec_X2b(abs(hnum(j))) = spec_X2b(abs(hnum(j))) + conj(spec_X2bZ(j));
        end
    end
end

dX2b = zeros(size(dX2));
for l = 1:length(w)
    dX2b = dX2b + real(AdX2b(l) * exp(-1i*w(l)*ta));
end

X2b = zeros(size(dX2b));
for l = 1:length(w)
    X2b = X2b + real(spec_X2b(l) * exp(-1i*w(l)*ta));
end

%%

X = X1 + X2 + X2b;

spec_X = spec_X1 + spec_X2 + spec_X2b;
phase_X = wrapToPi(atan2(real(spec_X), imag(spec_X)) - pi/2);

X_rct = zeros(size(X));
for j = 1:length(w)
    X_rct = X_rct + abs(spec_X(j)) .* cos(phase_X(j) + w(j) * ta);
end

%%

%{
fileName = 'displ_second_order.dat';
fileID = fopen(fileName,'w');
fprintf(fileID, '# time displacement\n');
for j = 1:length(ta)
    fprintf(fileID, '%15.10e %15.10e\n', ta(j), X(j));
end
fclose(fileID);
%}

%{
% Output 2nd order spectrum
% X = amp * cos(phase + omega * t)
fileName = 'spec_second_order.dat';
fileID = fopen(fileName,'w');
fprintf(fileID, '# amp omega phase\n');
for j = 1:length(w)
    fprintf(fileID, '%15.10e %15.10e %15.10e\n', abs(spec_X(j)), w(j), phase_X(j));
end
fclose(fileID);
%}

%%

fig1 = figure('pos',[100 100 1400 800]);

sub1 = subplot(3,2,1,'Parent',fig1);
grid on;
hold on;
plot(ta, Ea, '-b', 'Parent', sub1, 'LineWidth', 1);
plot(ta, Eb, '-g', 'Parent', sub1, 'LineWidth', 1);
hold off;
axis([-inf inf -inf inf]);
legend({'Ea' 'Eb'});
xlabel('t [s]');

sub2 = subplot(3,2,2,'Parent',fig1,'YScale','log');
grid on;
hold on;
plot(w, abs(Aa), '+-b', 'Parent', sub2);
plot(w, abs(Ab), 'x-g', 'Parent', sub2);
hold off;
axis([0 w(end) 1e-5 inf]);
set(gca,'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1]);
legend({'Aa' 'Ab'},'Location','NorthEastOutside');

sub3 = subplot(3,2,3,'Parent',fig1,'YScale','log');
grid on;
hold on;
plot(w, abs(AdX1), '+-b', 'Parent', sub3, 'LineWidth', 1);
plot(w, abs(AdX2), '-^c', 'Parent', sub3, 'LineWidth', 1);
plot(w, abs(AdX2b), 'x-g', 'Parent', sub3, 'LineWidth', 1);
hold off;
axis([0 w(end) 1e-5 inf]);
set(gca,'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1]);
legend({'spec dX1' 'spec dX2' 'spec dX2b'},'Location','NorthEastOutside');

sub4 = subplot(3,2,4,'Parent',fig1,'YScale','log');
grid on;
hold on;
plot(w, abs(spec_X1), '+-b', 'Parent', sub4, 'LineWidth', 1);
plot(w, abs(spec_X2), '^-c', 'Parent', sub4, 'LineWidth', 1);
plot(w, abs(spec_X2b), 'x-g', 'Parent', sub4, 'LineWidth', 1);
hold off;
axis([0 w(end) 1e-5 inf]);
set(gca,'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1]);
legend({'spec X1','spec X2','spec X2b'},'Location','NorthEastOutside');

sub5 = subplot(3,2,5,'Parent',fig1);
grid on;
hold on;
plot(ta, dX1, '-b', 'Parent', sub5, 'LineWidth', 1);
plot(ta, dX2, '-c', 'Parent', sub5, 'LineWidth', 1);
plot(ta, dX2b, '-g', 'Parent', sub5, 'LineWidth', 1);
plot(ta, dX1+dX2+dX2b, '-r', 'Parent', sub5, 'LineWidth', 1);
hold off;
axis([-inf inf -inf inf]);
legend({'dX1' 'dX2' 'dX2b' 'sum'});

sub6 = subplot(3,2,6,'Parent',fig1);
grid on;
hold on;
plot(ta, X, '-r', 'Parent', sub6, 'LineWidth', 1);
plot(ta, X_rct, '.k', 'Parent', sub6, 'LineWidth', 1);
hold off;
axis([-inf inf -inf inf]);
legend({'X' 'Xrct'});
