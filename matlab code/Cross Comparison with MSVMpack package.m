addpath('~/Documents/MSVMpack1.5/matlab/')
clear
load WW-3class.mat;

model = trainmsvm(X,Y, "-m WW -k 1 -c 1 -q -u -a 1");

model = loadmsvm( 'noname')


% Y
% [model.alpha'*X, model.b]


alpha = (sum(model.alpha')'*ones(1,4) - model.alpha).*Y_mat - model.alpha;
[alpha'*X, model.b]


X = [-2.0 1.0;
    2.0 1.0;
    -2.0 -1.0;
    -2.0 1.0;
    2.0 1.0;
    -2.0 -1.0;
    -2.0 1.0;
    2.0 1.0;
    -2.0 -1.0;
    -2.0 1.0;
    2.0 1.0;
    -2.0 -1.0]


Y = [1;2;3;1;2;3;1;2;3;1;2;3;1;2;3];
Y = double(Y);
Y_mat = [1 0 0;
          0 1 0;
          0 0 1;1 0 0;
          0 1 0;
          0 0 1;
          1 0 0;
          0 1 0;
          0 0 1;
          1 0 0;
          0 1 0;
          0 0 1]

model = trainmsvm(X,Y, "-m WW -k 1 -c 1 -u -a .9999");

model = loadmsvm( 'noname')

alpha = (sum(model.alpha')'*ones(1,3) - model.alpha).*Y_mat - model.alpha;
[alpha'*X, model.b]

