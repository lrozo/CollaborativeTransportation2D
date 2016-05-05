function [Mu, Sigma] = gaussianProduct(model, p)

% Leonel Rozo, 2014
% 
% Compute the product of Gaussians 

% GMM products 
for i = 1 : model.nbStates 
  SigmaTmp = zeros(model.nbVar);
  MuTmp = zeros(model.nbVar,1);
  
  for m = 1 : model.nbFrames 
    MuP = p(m).A * model.Mu(:,m,i) + p(m).b; 
    SigmaP = p(m).A * model.Sigma(:,:,m,i) * p(m).A'; 
    SigmaTmp = SigmaTmp + inv(SigmaP);
    MuTmp = MuTmp + SigmaP\MuP; 
  end
  Sigma(:,:,i) = inv(SigmaTmp);
  Mu(:,i) = Sigma(:,:,i) * MuTmp;    
end