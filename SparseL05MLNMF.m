function [Af, S] = SparseL05MLNMF(X, P, a0, tau, d, Lmax, Tmax)
% Multilayer NMF by Roozbeh Rajabi and Hassan Ghassemian
% ----------------------
% Email: rajabi@ieee.org
% ----------------------
% Reference: R. Rajabi, H. Ghassemian, "Spectral Unmixing of Hyperspectral
% Imagery Using Multilayer NMF," IEEE Geoscience and Remote Sensing Letters
% , vol.12, no.1, pp.38-42, Jan. 2015.
% ----------------------
% SparseL05MLNMF(X, P, a0, tau, d, Lmax, Tmax)
% Inputs:
% X (Observation): L*N 
% P (Number of signature)
% a0, tau (Regularization parameter)
% d (ASC)
% Lmax (Maximum number of layers)
% Tmax (Maximum number of iterations in each layer)
% Outputs:
% A (Signatures): L*P;
% S (Abundances): P*N;

[L,N] = size(X);

% VCA-FCLS Initialization
[A, ~] = VCA(X,'Endmembers', P,'verbose','on');

% ASC
Xb = [X; d*ones(1,N)];
Ab = [A; d*ones(1,P)];
S = (Ab'*Ab)^-1*Ab'*Xb;

Ap = ones(size(A));
count = 0;
Stemp = S;
Af = eye(L);
for l = 1:Lmax
    if l ~= 1
        A = rand(P,P);
    end
    for t = 1:Tmax
        aA = a0 * exp(-t/tau);
        aS = 2 * a0 * exp(-t/tau);
        % Updating Rules
        A = A.*(X*S')./(A*(S*S')+0.5*aA*A.^-0.5);
        A(A<eps) = eps;
        % ASC
        Xb = [X; d*ones(1,N)]; Ab = [A; d*ones(1,P)];
        S(S<eps) = eps;
        % S>=10^-4
        S1 = S.*(Ab'*Xb)./(Ab'*Ab*S+0.5*aS*S.^-0.5);
        % S<10^-4
        S2 = S.*(Ab'*Xb)./(Ab'*Ab*S);
        Stemp(S>=10^-4) = S1(S>=10^-4);
        Stemp(S<10^-4) = S2(S<10^-4);
        S = Stemp;
        if t ~= 1
            err = norm(A-Ap);
            disp([num2str(err) '@' num2str(t)])
            if(err < 1.0e-4)
                count = count + 1;
                if (count == 10)
                    disp(['done!' 'Tolerance:' num2str(err) '@' num2str(t)]);
                    break;
                end
            else
                count = 0;
            end
        end
        Ap = A;
    end
    if t == Tmax
        disp(['done!' 'Tolerance:' num2str(err) '@' num2str(t)]);
    end
    X = S;
    Af = Af * A;
end
