clc; clear all; close all;

% Load Synthetic Data
load syntheticData;
SNR = 25;

% Add Noise
variance = sum(mixed(:).^2)/10^(SNR/10)/M/N/D;
n = sqrt(variance)*randn([D M*N]);
mixed = mixed' + n;
clear n;

% remove noise
[UU, SS, VV] = svds(mixed,P);
Lowmixed = UU'*mixed;
mixed = UU*Lowmixed;
EM = UU'*Ao;

X = mixed;

% Parameters
a0 = 0.1;
tau = 25;
d = 25;
Lmax = 10; Tmax = 400;

% MLNMF Unmixing
tic
[Ae, Se] = SparseL05MLNMF(X, P, a0, tau, d, Lmax, Tmax);
toc

% Permutation
perm = permute_corr(Ao,Ae);
Ae = Ae * perm;
% Rescaling
Ae = Ae./repmat(max(Ae), size(Ae,1), 1);
Ae = Ae.*repmat(max(Ao), size(Ae,1), 1);

% Plot signatures
figure; hold on;

plot(wlen,Ao(:,1:P),'LineWidth',2.25)
plot(wlen,Ae(:,1:P),'-.','LineWidth',1.25)

xlabel('Wavelenghts(\mum)')
ylabel('Reflectance')
legend(name)
axis([wlen(1) wlen(end) 0 1]);

% SAD & AAD
[sad, rmsSAD] = AD(Ao',Ae');
disp(sad')
disp(mean(sad));

Se = perm' * Se;
[aad, rmsAAD] = AD(abf',Se');
disp(rmsAAD);
