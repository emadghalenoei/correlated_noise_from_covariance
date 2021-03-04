% How to make correlated noise to data from a covariance matrix

% more info can be found at my journal paper entitled "Gravity and Magnetic Joint Inversion to 
                                   ... Resolve Basement and Salt Structures with Reversible-Jump Algorithm"
                                       
% descibing inputs, processing and outputs
% inputs: x and data
% outputs: noisy data

% Emad Ghalenoei, Univesity of Calgary, March 2021
%
% See also TheOtherFunction.
% Copyright: This is published under GJI manuscript entitled "Gravity and Magnetic Joint Inversion to 
                                   ... Resolve Basement and Salt Structures with Reversible-Jump Algorithm"
% You are allowed to use, as long as you cite this work.


clc
clear
close all

x =  (0:0.1:2*pi)';
data = sin(x);

noise_level = 0.1; % 10%
sigma = noise_level.*max(abs(data)); % noise standard deviation 

ind=0:1:length(data)-1;
Dij=toeplitz(ind);
% create a sinusoidal correlation matrix (R)
R = exp(-Dij./20).*cos(2*pi.*Dij./20);
Cov = (sigma^2) * R; % cov matrix

np = normrnd(0,1,size(data)); % generate random white noise N(0,1)
L = chol(Cov,'lower'); % lower cholesky matrix
correlated_noise = L * np;
noisy_data = data + correlated_noise;

[autocor,lags] = xcov(correlated_noise,'coeff');
autocor = autocor(lags>=0);

figure
subplot(2,1,1)
plot(x,data,'LineWidth',2)
hold on
plot(x,noisy_data,'k.-','LineWidth',2)
xlabel('x'); ylabel('data'); legend('data','noisy data')
set(gca,'fontsize',15,'fontweight','bold')

subplot(2,1,2)
plot(lags(lags>=0),autocor,'LineWidth',2)
xlabel('x'); ylabel('auto correlation of noise');
set(gca,'fontsize',15,'fontweight','bold')