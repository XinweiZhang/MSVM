addpath('~/Documents/MSVMpack1.5/matlab/')
clear
load ../WW-4class.mat;

model = trainmsvm(X,Y, "-m WW -k 1 -c 1 -q");

model = loadmsvm( 'noname')


Y
[model.alpha'*X, model.b]


alpha = (sum(model.alpha')'*ones(1,4) - model.alpha).*Y_mat - model.alpha
[alpha'*X, model.b]

