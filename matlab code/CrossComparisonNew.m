addpath('~/Documents/MSVMpack1.5/matlab/')
clear
load data.mat;

% model = trainmsvm(X,y, "-m CS -k 1 -c 1 -u ");
% Y
% alpha = sum((model.alpha)')'*ones(1,6).*Y - model.alpha
% 
% [model.alpha'*X, model.b]
% 
% [alpha'*X, model.b]


% RBF kernel
model = trainmsvm(X,y, "-m CS -k 1 -p 1 -c 1 -q -u -a .98");

model = loadmsvm( 'noname')
(model.alpha'*X)'

% model = trainmsvm(X,Y, "-m WW -k 1 -c 1 -q -u -a 1");

model.alpha(1:5,:)
%  

% % % 
alpha = (sum(model.alpha')'*ones(1,m) - model.alpha).*Y - model.alpha;

alpha(1:5,:)

[alpha'*K, model.b]'

% 
% 
% X = [-2.0 1.0;
%     2.0 1.0;
%     -2.0 -1.0;
%     -2.0 1.0;
%     2.0 1.0;
%     -2.0 -1.0;
%     -2.0 1.0;
%     2.0 1.0;
%     -2.0 -1.0;
%     -2.0 1.0;
%     2.0 1.0;
%     -2.0 -1.0]
% 
% 
% Y = [1;2;3;1;2;3;1;2;3;1;2;3;1;2;3];
% Y = double(Y);
% Y_mat = [1 0 0;
%           0 1 0;
%           0 0 1;1 0 0;
%           0 1 0;
%           0 0 1;
%           1 0 0;
%           0 1 0;
%           0 0 1;
%           1 0 0;
%           0 1 0;
%           0 0 1]
% 
% model = trainmsvm(X,Y, "-m WW -k 1 -c 1 -u -a .9999");
% 
% model = loadmsvm( 'noname')
% 
% alpha = (sum(model.alpha')'*ones(1,3) - model.alpha).*Y_mat - model.alpha;
% [alpha'*X, model.b]
% 