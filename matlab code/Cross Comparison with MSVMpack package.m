addpath('~/Documents/MSVMpack1.5/matlab/')
clear
load WW-4class.mat;

model = trainmsvm(X,Y, "-m WW -k 1 -c 1 -q -u");

model = loadmsvm( 'noname')


% Y
% [model.alpha'*X, model.b]


alpha = (sum(model.alpha')'*ones(1,4) - model.alpha).*Y_mat - model.alpha;
[alpha'*X, model.b]


X = [-2 1;
    2 0.9;
    2 1;
    2 0.9;
    -2 -1]
X = double(X);
Y = [1;1;2;2;3];
Y = double(Y);
Y_mat = [1 0 0;
          0 1 0;
          0 0 1]

model = trainmsvm(X,Y, "-m WW -k 1 -c 100 -u");

model = loadmsvm( 'noname')

alpha = (sum(model.alpha')'*ones(1,4) - model.alpha).*Y_mat - model.alpha;
[alpha'*X, model.b]



