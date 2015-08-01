function [A, rho, Pz_x, zRange] = em_new_indivzrange(x,A,rho,kGam,thGam,sigmaX,varargin)

% Usage: [A, rho, Pz_x] = em_new_indivzrange(x,A,rho,kGam,thGam,sigmaX,varargin)
%
% The function performs an EM step on the data batch x. The scalar-vauled z
% comes from a Gamma prior (w parameters thGam & kGam), while latent u
% comes from a Normal prior (w mean 0 & covariance rho). Observation noise
% is sigmaX.
%
% The function is the same as em_new but intends to set the range of z values
% automatically instead of sweeping through a whole range of values for z's
%
% Input arguments:
% x     data; dimensions: Dim(x) x Batch size
% A     feature matrix; dimensions: Dim(x) x Dim(u)
% rho   the covariance matrix of the normally-distributed latent, it is
%       referred to as C in the derivation; dimensions: Dim(u) x Dim(u)
% varargin sets the value of rhoOnly: if set to 1, then learning of A is
%       omitted and a nan is given for the output A
%
% Output arguments: 
% A     feature matrix
% rho   the covariance matrix of the normally-distributed latent, it is
%       referred to as C in the derivation
% Pz_x  posterior of the scalar-valuad z for each stimulus. Maximum
%       z/number of points evaluated might change depending on the
%       posterior mean's location. This effectively implies that the number
%       of mixture components might change
% zRange the values of z where the Pz_x has been evaluated
%
% Written by Gergo Orban (go223@cam.ac.uk)

args=varargin;
nargs=length(args);
if (nargs>0)
    rhoOnly=args{1};
else
    rhoOnly=0;
end
if (nargs>1),
    Dz = args{2};
else
    Dz = 100;
end

[Dx T]=size(x);
Du=size(rho,1);

rhoInv = inv(rho);
ATA = A' * A;
ACAT = A * rho * A';
Sigma0 = eye(Dx)*(sigmaX^2);
x0=zeros(Dx,1);

zRange = zeros(1,T, Dz);
zPostMaxs=zeros(1,T);
zRanges=zeros(2,T);
for tt=1:T
    [zPostMaxs(tt) zRanges(1,tt) zRanges(2,tt)]= ...
        get_pz_x_max(x(:,tt), Sigma0, ACAT, x0,kGam,thGam);
    if (Dz>1)
        zRange(1,tt,:) = linspace(zRanges(1,tt), zRanges(2,tt),Dz);
    else
        zRange(1,tt,1) = zPostMaxs(tt); 
    end
end;

lPz_x = ones(1,T,Dz)*(-1000000000000);
Pz_x = zeros(1,T,Dz);

Ksi_u=zeros(Du,Du);
Ksi_y=zeros(Du,Du);
Psi_y=zeros(Dx,Du);

for tt=1:T,
    if (Dz>1),
        for zi=1:Dz,
            lPz_x(:,tt,zi) = log(gampdf(zRange(1,tt,zi),kGam,thGam)) + ...
                normpdfln(x(:,tt), x0, [], Sigma0+zRange(1,tt,zi)^2 * ACAT);     % 1*T*z
        end;
        Pz_x(1,tt,:) = exp(lPz_x(1,tt,:) - max(lPz_x(1,tt,:)));
        Pz_x(1,tt,:) = Pz_x(1,tt,:) / sum(Pz_x(1,tt,:));
    else
        Pz_x(1,tt,1)=1;
        lPz_x(1,tt,1)=0;
    end
    
    Psi_y_x = zeros(1,Du);
    for zi=1:Dz,
        Sigma = inv(rhoInv + zRange(1,tt,zi)^2/(sigmaX^2) * ATA);   % Du * Du
        mu_x = zRange(1,tt,zi) /(sigmaX^2) * Sigma * A';        % Du * Dx
        mu = mu_x * x(:,tt);                                    % Du * 1
        Ksi_uact = Pz_x(1,tt,zi) * ((mu * mu') + Sigma);
        Ksi_u = Ksi_u + Ksi_uact ;
                
        if (rhoOnly == 0)
            Ksi_y = Ksi_y + zRange(1,tt,zi)^2 * Ksi_uact;
            Psi_y_x = Psi_y_x + Pz_x(1,tt,zi) * zRange(1,tt,zi) * mu';
        end
    end;
    if (rhoOnly == 0),
        Psi_y = Psi_y + x(:,tt) * Psi_y_x;
    end
end
rho = Ksi_u / T;

if (rhoOnly==0),
    A = Psi_y / Ksi_y;
else
    A=NaN;
end
