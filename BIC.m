function [BIC1, penalty] = BIC(data, model, L, N)

% Leonel Rozo, 2015
% 
% This function computes the Bayesian Information Criterion (BIC) score for
% a given model \lambda={\pi, \miu, \Sigma}. The model can be a TP-GMM or
% some similar approaches (e.g., HMM). It outputs the BIC scores following 
% the Calinon formula.
%
% input:
%   matrix      data    Training data
%   structure   model   GMM model (priors, means and covariances)
%   float       L       Likelihood previously computed

% Dimension of datapoints
D = size(data, 1 ) ;
% Number of gaussian components K
K = model.nbStates ;
% Number of frames
P = model.nbFrames;


%% * * * * * * * * * * Bayesian Information Criterion * * * * * * * * * * *

% Computing number of free parameters required for a GMM of K components
% assuming full covariance matrices
Np = (K - 1) + K * P * (D + 0.5 * D * (D+1)) ;
    
% Computing BIC's second term - Penalty factor
penalty = 0.5 * (Np) * log(N);

% Computing BIC_1
BIC1 = -L + penalty ;
