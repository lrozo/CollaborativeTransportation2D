function pHRI_toyExample_PlanarManip01

% Leonel Rozo, 2015
%
% This code shows a simple human-robot cooperative transportation task in a
% planar scenario. The task consists of transporting an object from an 
% initial position to a target location.
%
% This code shows:
% 
% 1. Computation of virtual attractor trajectories from the dynamics
%    observed during the demonstrations of the collaborative task.
% 2. TP-GMM learning of the collaborative task.
% 3. BIC-based model selection.
% 4. Stiffness estimation built on a convex optimization.
% 5. Reproduction using GMR with adaption to new configurations of the task
%    parameters.
%
% * TASK PARAMETERS
%   Three different task parameters are defined in this experiment,
%   namely: 
%   - Initial position and orientation of the object
%   - Target position and orientation of the object
%   - Irrelevant frame that randonmly varies its position and orientation
%
% * DEPENDENCIES:
%   - CVX   Matlab toolbox for convex optimization.(http://cvxr.com/cvx/)

%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -----> Model variables
model.nbStates = 4;    % Number of Gaussians in the model
model.nbFrames = 3;    % Number of candidate frames of reference
model.nbVar    = 3;    % Dimension of the datapoints [t y1 y2]
model.kP       = 600;  % Initial stiffness gain 
model.kV       = 50;   % Damping gain 
model.dt       = 0.01; % Time step

% -----> Variables for controlling different parts of the code
needsModel   = 0;   % For learning or loading a TP-GMM 
BICselection = 0;   % For carrying out a BIC-based model selection
stiffEstimate= 0;   % For estimating/loading stiffness matrices
saveModelMat = 1;   % Save TP-GMM for MATLAB
saveStiffEst = 1;   % Save stiffness estimations

% ------ Files and folders
dataPath  = 'data/'; % Data path


%% Load demonstrations and task parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([dataPath 'Data.mat']);
nbData    = size(s(1).DataP,2); % Number of datapoints
demons    = 1:5;                % Demonstrations to train the model
nbSamples = length(demons);     % Number of demonstrations
% Ids
posId = (1:2); velId = (3:4); accId = (5:6);

% ------ CREATING TRAINING DATA [t y1 y2]
for n = demons
  s(n).nbData = nbData; 
  s(n).Data   = zeros(model.nbVar, s(n).nbData);   
  % Attractor's path Eq. (3)
  % y = \ddot{x}*(Kp^{-1}) + \dot{x}*(Kp^{-1}*Kv) - Kp^{-1}*Fs + x 
  s(n).Data = [ [1:nbData]*model.dt; ...
    (s(n).DataP(accId,:)*(1/model.kP) + ...    
    s(n).DataP(velId,:)*(model.kV/model.kP) - s(n).DataF*0.0 +...
    s(n).DataP(posId,:)) ];
end

% ------ CREATING TENSOR
% Project data to every frame and create the tensor
Data=[];
cnt = 1;
for n = demons  
  for m = 1 : model.nbFrames
    % [ t y1 y2 ] projected at each frame m  
    tt1 = 1+nbData*(cnt-1);
    tt2 = cnt*nbData;
    Data(:,m,tt1:tt2) = s(n).p(m).A \ ...
      (s(n).Data - repmat(s(n).p(m).b,1,s(n).nbData));     
  end  
  cnt = cnt+1;
end


%% Learning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Parameters estimation of TP-GMM:\n');
if(needsModel)
  % Train several models and select the one with the lowest BIC  
  if(BICselection) 
    fprintf('BIC-based model selection...\n');    
    K = 1:10;       % Vector of number of states to try
    BICvalues = []; % BIC value for each model
    penalty   = []; % "Penalty" values computed at BIC 
    LogL      = []; % Log-likelihood
    m         = []; 
    
    for k = K  
      model.nbStates = k;
      model = init_tensorGMM_timeBased(Data, model); %Initialization 
      [model, ~] = EM_tensorGMM(Data, model); % Training Eq. (5)-(8)     
      m(k).model = model;
      
      % The method computes the log-likelihood for every GMM obtained
      % for every set of task-parameters of the demonstrations. The TP-GMM
      % was previously trained.
      for n = demons        
        % Computing resulting GMM  
        [Mu, Sigma] = gaussianProduct(model, s(n).p);
  
        for i = 1 : model.nbStates
          % Computing likelihood
          Lklhd(i,:) = model.Priors(i) * gaussPDF(s(n).Data, Mu(:,i), ...
            Sigma(:, :, i));          
        end
        %Compute log-likelihood
        LL(n) = sum(log(sum(Lklhd,1)));
        % BIC for current demo
        [auxBIC(n), auxPenalty(n)] = BIC(s(n).Data, model, LL(n), nbData);
      end
      % Storing mean values
      LogL = [LogL mean(LL)];
      BICvalues = [BICvalues mean(auxBIC)];
      penalty = [penalty mean(auxPenalty)];
    end
    % Choosing the model with the minimum BIC
    [~, minInd] = min(BICvalues);
    model = m(K(minInd)).model; 
    fprintf('Done!\n');    
  else  
    model = init_tensorGMM_timeBased(Data, model); %Initialization   
    model = EM_tensorGMM(Data, model);
  end
  
  if(saveModelMat)
    save([dataPath 'model.mat'], 'model');
  end
else
  fprintf('Loading model...');  
  load([dataPath 'model.mat']);
  fprintf('Done!\n');  
end


%% Stiffness estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(stiffEstimate)
  fprintf('Stiffness matrices estimation...');  
  Y = []; W = []; X = [];
  for n = demons
    w = [];
    auxX = [];    
  
    % Computing resulting GMM  
    [Mu, Sigma] = gaussianProduct(model, s(n).p);
  
    for i = 1 : model.nbStates
      % Computing influence weights Eq. (5)
      w(:,i) = gaussPDF(s(n).Data(1,:), Mu(1,i), Sigma(1, 1, i));
      % Computing x tilde 
      auxX(:,:,i) = repmat(Mu(2:3,i), 1, nbData) - s(n).DataP(posId,:);      
    end  
    % Weights normalization
    auxW = w ./ repmat(sum(w, 2), 1, model.nbStates);
    % Computing variable v 
    auxY = s(n).DataP(accId,:) + model.kV * s(n).DataP(velId,:) - s(n).DataF;    
    % Concatenating variables
    W = [W ; auxW]; X = [X auxX]; Y = [Y auxY];  
  end

  for i = 1 : model.nbStates
    % Weighting data Eq. (12) 
    XtW = (X(:,:,i)' .* repmat(W(:, i), 1, size(X,1)))';
    YtW = Y .* repmat(W(:, i)', size(Y,1), 1);
     
    % ------ Convex optimization Eq. (14)
    d = length(posId);
    cvx_begin sdp
      variable Kcvx(d,d) symmetric;
      minimize( norm(Kcvx * XtW - YtW, 2)); % Eq. (13)
      Kcvx >= eye(d);
    cvx_end
    KpCVX(:,:,i) = Kcvx
  end
  
  if( saveStiffEst )
    save([dataPath 'KpCVX.mat'], 'KpCVX');      
  end
  fprintf('Done!\n');  
else
  fprintf('Loading stiffness matrices...');  
  load([dataPath 'KpCVX.mat']);
  fprintf('Done!\n');  
end
% ------Clamping the estimated Kp for each state
% Rescale Kp to stay within a specific range [kPmin, kPmax]. Useful for
% practical implementations on real robots
KpMax = 1800;
KpMin = 400;
for k = 1 : model.nbStates
  % Using only the diagonal terms  
  Kp(:,:,k) = diag(diag(KpCVX(:,:,k)));
  % Eigencomponents decomposition
  [V(:,:,k), Dtmp1] = eig(Kp(:,:,k));   
  lambda(:,k) = diag(Dtmp1);   
end

lambda_min = min(min(lambda));
lambda_max = max(max(lambda));

for k = 1 : model.nbStates
  % Rescale each eigenvalue such that they lie in the range [kPmin,kPmax]
  lambdaFactor = (lambda(:,k) - lambda_min) ./ (lambda_max-lambda_min);
  Dtmp1 = diag(((KpMax-KpMin) .* lambdaFactor) + KpMin);
  % Reconstruction from the modified eigencomponents
  Kp(:,:,k) = (V(:,:,k) * Dtmp1) / V(:,:,k);       
end


%% Reproduction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
repros = [6:8]; % Data2
in     = 1;     % Input variables
out    = [2 3]; % Output variables

model.estKp = Kp;

for n = repros
  % Setting initial Cartesian position of the robot  
  x0 = s(n).DataP(posId,1);
  input = [1:nbData]*model.dt;
  
  % DS-GMR
  r(n) = reproduction_tensorGMM_DS_Kp(input, model, s(n), x0, in, out,...
    s(n).DataF);    
end


%% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting BIC-based model selection results
if(needsModel && BICselection)
  figure('position',[20,50,1100,600], 'color', [1 1 1]);
  subplot(2,2,1); hold on; box on;   % BIC values
    plot(BICvalues, '-b', 'LineWidth', 2);
    plot(BICvalues, 'or', 'MarkerSize', 8, 'MarkerFaceColor', 'r' ) ;
    set(gca, 'XTick', 1:length(K), 'XTickLabel', K);
    set(gca, 'FontSize', 22);
    xlabel('$K$', 'Fontsize', 18,'interpreter','latex');
    ylabel('BIC', 'Fontsize', 18,'interpreter','latex');
end
      
% Demonstrations, local models, and reproductions
figure('PaperPosition', [0 0 18 4.5], 'position', [100,100,1400,900],...
  'color',[1 1 1]); 
% Colors
GMMclrs = [[.95 .10 .10] ; [.12 .12 .99] ; [.95 .56 .12] ; [0.15 0.95 0.45]];
% Labels
xLbls = ['$x_1$' ; '$x_2$'];
fLbls = ['$f_1$' ; '$f_2$'];
ySLbls = ['$y_1^{\mathcal{S}}$';'$y_2^{\mathcal{S}}$'];
yTLbls = ['$y_1^{\mathcal{T}}$';'$y_2^{\mathcal{T}}$'];
yILbls = ['$y_1^{\mathcal{I}}$';'$y_2^{\mathcal{I}}$'];
frmLbls(1).Lbls = ySLbls; frmLbls(2).Lbls = yTLbls; frmLbls(3).Lbls = yILbls;
% Axes limits
xLims(:,:,1) = [[-0.05 2.0] ; [-0.50 2.0]];
xLims(:,:,2) = [[-2.0 0.5] ; [-2.0 0.5]];
xLims(:,:,3) = [[-1.0 1.5] ; [-1.5 1.0]];

% Local models and demonstrations
subplot(2,4,1); hold on; box on;
  for n = demons
    for m = 1 : model.nbFrames  
      plot2Dframe(s(n).p(m,1).A(2,2:3), s(n).p(m,1).A(3,2:3),...
        s(n).p(m,1).b(2:3), GMMclrs(m,:), GMMclrs(m,:), 0.2);  
    end  
    plot(s(n).DataP(1,:), s(n).DataP(2,:), '-', 'color', [0.3 0.3 0.3]); % Robot position
  end
  xlabel('$x_1$', 'Fontsize', 18, 'interpreter', 'latex'); 
    ylabel('$x_2$', 'Fontsize', 18, 'interpreter', 'latex'); 
    axis equal; axis square; 
    set(gca, 'Fontsize', 16);
    set(gca,'Layer','top');
    axis([-1.0 1.0 -1.0 1.0]);
text(-0.7,1.25,'Demonstrations','Fontsize', 18, 'interpreter', 'latex'); 

for m = 1 : model.nbFrames
  cnt = 1;  
  for n = 1 : nbSamples
    % Auxiliar indexes
    tt1 = 1+nbData*(cnt-1);
    tt2 = cnt*nbData;  
    subplot(2,4,m+1); hold on; box on;     
    if n == 1
      plotGMM(reshape(model.Mu(2:3,m,:),[2 model.nbStates]),...
        reshape(model.Sigma(2:3,2:3,m,:), [2 2 model.nbStates]),...
        GMMclrs(m,:));
      plot2Dframe([1 0], [0 1], zeros(1,2), GMMclrs(m,:), GMMclrs(m,:), 0.2);          
    end
    plot(reshape(Data(2,m,tt1:tt2), [1 nbData]), ... 
      reshape(Data(3,m,tt1:tt2), [1 nbData]), '-', ...
      'color', [0.3 0.3 0.3]); % Robot attractor
    plot(Data(2,m,tt1), Data(3,m,tt1), '.', 'Markersize', 18, ...
      'color', [0.0 0.0 0.0]); % Start pos
    plot(Data(2,m,tt2), Data(3,m,tt2), 'x', 'Markersize', 10, ...
      'color', [0.0 0.0 0.0]); % Target pos
    xlabel(frmLbls(m).Lbls(1,:), 'Fontsize', 18, 'interpreter', 'latex'); 
    ylabel(frmLbls(m).Lbls(2,:), 'Fontsize', 18, 'interpreter', 'latex'); 
    axis equal; axis square; 
    set(gca, 'Fontsize', 16);
    set(gca,'Layer','top');
    axis([xLims(1,:,m) xLims(2,:,m)]);
    cnt = cnt+1;
  end
end
text(-5.2,1.3,'Attractor projected on the task parameters',...
    'Fontsize', 18, 'interpreter', 'latex'); 

% Resulting model and reproductions
cnt = 1;
for n = repros
  subplot(2,3,cnt+3); hold on; box on;     
    plotGMM(r(n).Mu(2:3,:), r(n).Sigma(2:3,2:3,:), GMMclrs(4,:));
    plot2Dframe([1 0], [0 1], zeros(1,2), GMMclrs(4,:), GMMclrs(4,:), 0.2);          
    plot2Dframe([1 0], [0 1], s(n).p(1).b(2:3), GMMclrs(1,:), ...
      GMMclrs(1,:), 0.2);
    plot2Dframe([1 0], [0 1], s(n).p(2).b(2:3), GMMclrs(2,:), ...
      GMMclrs(2,:), 0.2);
    % Robot attractor
    h1 = plot(r(n).Data(8,:), r(n).Data(9,:), '--', 'color', [0.7 0.7 0.7]); 
    % Plotting sensed force
    for t = 1 : nbData        
      h2 = quiver(r(n).Data(2,t), r(n).Data(3,t), s(n).DataF(1,t), ...
        s(n).DataF(2,t), 0.02, 'color', [0.8 0.8 0.6]);          
    end    
    % Robot position
    h3 = plot(r(n).Data(2,:), r(n).Data(3,:), '-', 'color', [0.3 0.3 0.3]); 
    plot(r(n).Data(2,1), r(n).Data(3,1), '.', 'Markersize', 18, ...
      'color', [0.0 0.0 0.0]); % Start pos
    plot(r(n).Data(2,end), r(n).Data(3,end), 'x', 'Markersize', 10, ...
      'color', [0.0 0.0 0.0]); % Target pos
    xlabel('$x_1$', 'Fontsize', 18, 'interpreter', 'latex'); 
    ylabel('$x_2$', 'Fontsize', 18, 'interpreter', 'latex'); 
    axis equal; axis square; 
    set(gca, 'Fontsize', 16);
    set(gca,'Layer','top');
    axis([-1.1 1.1 -1.1 1.1]);
    
    cnt = cnt + 1;
end
text(-4.2,1.2,'Reproductions for different task parameters',...
    'Fontsize', 18, 'interpreter', 'latex'); 
legend([h1, h2, h3], 'Attractor', 'Sensed force', 'Robot pos', ...
  'Location','southeast');

% Reproduction data, activation of Gaussian components, and stiffness
figure('PaperPosition', [0 0 18 4.5], 'position', [100,100,1500,600],...
  'color',[1 1 1]);
xx = round(linspace(1, 64, model.nbStates));
clrmap = colormap('Jet');
clrmap = min(clrmap(xx,:), .95);

for i = 1 : 2
  subplot(2,3,i); hold on; box on; % Robot position and attractor          
    plot(r(n).Data(1,:), r(n).Data(7+i,:), '--', 'color', [0.7 0.7 0.7]); % Robot attractor
    plot(r(n).Data(1,:), r(n).Data(1+i,:), '-', 'color', [0.3 0.3 0.3]); % Robot position    
    ylabel(xLbls(i,:), 'Fontsize', 18, 'interpreter', 'latex'); 
    set(gca, 'Fontsize', 16);
    axis([0.0 1.0 -1.0 1.0]);
    if(i==1)
      legend('Attractor','Robot pos','Location','northwest');
    end
    
  subplot(2,3,i+3); hold on; box on; % Sensed forces
    plot(r(n).Data(1,:), r(n).Data(9+i,:), '-', 'color', [0.3 0.3 0.3]); 
    xlabel('$t$', 'Fontsize', 18, 'interpreter', 'latex'); 
    ylabel(fLbls(i,:), 'Fontsize', 18, 'interpreter', 'latex'); 
    set(gca, 'Fontsize', 16);
    axis([0.0 1.0 -7.0 7.0]);    
end
text(-0.3,28,'Evolution over time of one reproduction of the task',...
    'Fontsize', 18, 'interpreter', 'latex'); 
subplot(2,3,3); hold on; box on;       
  for i = 1 : model.nbStates
    plot(input, r(repros(1)).H(i,:), 'color', clrmap(i,:), 'LineWidth', 4);   
  end
  axis([0.02 1 -0.15 1.15]);
  ylim([-0.15 1.15]);     
  ylabel('$\gamma$', 'Fontsize', 20, 'interpreter', 'latex');         
  set(gca, 'Fontsize', 16);                  

subplot(2,3,6); hold on; box on;         
  plot(input, reshape(r(repros(1)).Kp(1,1,:), [1 size(r(repros(1)).Kp,3)]), ...
    'r-', 'LineWidth', 4);  
  plot(input, reshape(r(repros(1)).Kp(2,2,:), [1 size(r(repros(1)).Kp,3)]), ...
    'g--', 'LineWidth', 4);    
  axis([0.02 1 KpMin-100 KpMax+100]);
  xlabel('$t$', 'Fontsize', 20, 'interpreter', 'latex'); 
  ylabel('$K^{\mathcal{P}}$', 'Fontsize', 20, 'interpreter', 'latex');         
  set(gca, 'Fontsize', 16);
  legend('x_1','x_2','Location','northwest');
  