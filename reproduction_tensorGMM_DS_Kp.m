function r = reproduction_tensorGMM_DS_Kp(DataIn, model, rr, currPos, ...
             in, out, DataF_r)

% Leonel Rozo, 2014
%
% This function implements a DS-GMR using a stiffness matrix computed from
% the stiffnesses association to every Gaussian component of a GMM. This
% can be used when a stiffness estimation process has been carried out
% beforehand, so that every state of the model has a stiffness ellipsoid.
%
% Note that the code assumes that the task parameters are fixed along the
% repro. 
%

nbData = size(DataIn,2);

r.p = rr.p;

%------> GMM products (see Eq. (3) Calinon et al, Humanoids'2012 paper)
for i = 1 : model.nbStates 
  SigmaTmp = zeros(model.nbVar);
  MuTmp = zeros(model.nbVar,1);
  for m = 1 : model.nbFrames 
    MuP = r.p(m).A * model.Mu(:,m,i) + r.p(m).b; 
    SigmaP = r.p(m).A * model.Sigma(:,:,m,i) * r.p(m).A'; 
    SigmaTmp = SigmaTmp + inv(SigmaP);
    MuTmp = MuTmp + SigmaP\MuP; 
  end
  r.Sigma(:,:,i) = inv(SigmaTmp);
  r.Mu(:,i) = r.Sigma(:,:,i) * MuTmp;
end

%------> GMR
currVel = zeros(length(currPos), 1);
for t = 1 : nbData
  % Compute activation weight
  for i = 1 : model.nbStates
    r.H(i,t) = model.Priors(i) * gaussPDF(DataIn(:,t), r.Mu(in,i),...
      r.Sigma(in,in,i)); 
  end
  r.H(:,t) = r.H(:,t)/sum(r.H(:,t));
  
  % Evaluate the current target (Eq. (2) Calinon et al, Humanoids'2012),
  % and the stiffness of the virtual DS
  tarTmp = zeros(length(out),1);
  Kp = zeros(length(out));
  for i = 1 : model.nbStates
    % Attractor  
    tarTmp = tarTmp + r.H(i,t) * (r.Mu(out,i) + ...
      r.Sigma(out,in,i)/(r.Sigma(in,in,i)) * (DataIn(:,t)-r.Mu(in,i))); 
    % Stiffness
    Kp = Kp + r.H(i,t) * model.estKp(:,:,i);
  end
  r.Kp(:,:,t) = Kp;
  
  if(nargin == 7)
    currAcc = r.Kp(:,:,t)*(tarTmp-currPos) - model.kV*currVel + DataF_r(:,t);
  else
    currAcc = r.Kp(:,:,t)*(tarTmp-currPos) - model.kV*currVel;  
  end
  currVel = currVel + currAcc*model.dt;
  currPos = currPos + currVel*model.dt;
  
  % Keep a trace of data
  if(nargin == 7)
    r.Data(:,t) = [DataIn(:,t); currPos; currVel; currAcc; tarTmp; DataF_r(:,t)];
  else
    r.Data(:,t) = [DataIn(:,t); currPos; currVel; currAcc; tarTmp];  
  end
end