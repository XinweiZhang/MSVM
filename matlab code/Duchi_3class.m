clear
load ../Duchi-3class.mat;

n1 = size(X1);
n1 = n1(1);
n2 = size(X2);
n2 = n2(1);
n3 = size(X3);
n3 = n3(1);

cvx_begin quiet
    variables w1(p) w2(p) w3(p) b1 b2 b3;
    variables slack1(n1,m) slack2(n2,m) slack3(n3,m);
    variables xi1(n1) xi2(n2) xi3(n3); 

    minimize( sum_square([w1;w2;w3])/2 + C* sum([xi1;xi2;xi3]))
    subject to
        w1 + w2 + w3 == 0;
        b1 + b2 + b3 == 0;
        slack1(1:n1,1) == 1;
        X1*(w1 - w2) + b1 - b2 >= 1-slack1(:,2);
        X1*(w1 - w3) + b1 - b3 >= 1-slack1(:,3);
        slack2(:,2) == 1;
        X2*(w2 - w1) + b2 - b1 >= 1-slack2(:,1);
        X2*(w2 - w3) + b2 - b3 >= 1-slack2(:,3);
        slack3(:,3) == 1;
        X3*(w3 - w1) + b3 - b1 >= 1-slack3(:,1);
        X3*(w3 - w2) + b3 - b2 >= 1-slack3(:,2);
        slack1 >= 0;
        slack2 >= 0;
        slack3 >= 0;
        xi1 >= norms_largest(slack1,1,2)-1; 
        xi2 >= norms_largest(slack2,1,2)-1;
        xi3 >= norms_largest(slack3,1,2)-1;
        xi1 >= norms_largest(slack1,2,2)/2-1/2;
        xi2 >= norms_largest(slack2,2,2)/2-1/2;
        xi3 >= norms_largest(slack3,2,2)/2-1/2;
        xi1 >= norms_largest(slack1,3,2)/3-1/3;
        xi2 >= norms_largest(slack2,3,2)/3-1/3;
        xi3 >= norms_largest(slack3,3,2)/3-1/3;
cvx_end

[w1.' b1; w2.' b2; w3.' b3]


n = n1 + n2 + n3;


%%%%%%%%%%%%%%%%%4
% a = (rand(n,m)); % your matrix
% b = repmat(a,[1,1,m]);
% permute(b,[1,3,2])
% b(:,:,1)
% b(:,:,2)
cvx_begin
    variables w(m,p) b(m) slack(n,m) xi(n) t(n,m) u(n,m,m);
    minimize( sum(sum_square(w))/2 + C* sum(xi))
    subject to
        sum(w) == 0;
        sum(b) == 0;     
        ([X, ones(n,1)]*[w,b]'.*Y_mat)*ones(m,m) - [X, ones(n,1)]*[w,b]' >= (1-Y_mat).*(1-slack);
        slack >= 0;
        sum(slack.*Y_mat,2)==1;
        xi*ones(1,m) >= t + (ones(n,1)*1./(1:m)).*sum(u,3) - (ones(n,1)*1./(1:m)),[],dim=2;
%         xi>= max(t + (ones(n,1)*1./(1:m)).*sum(u,3)  - (ones(n,1)*1./(1:m)),[],2);
        repmat(t,[1,1,m]) + u >= permute(repmat(slack,[1,1,m]),[1,3,2]);
        u >= 0;
cvx_end



[w1.' b1; w2.' b2; w3.' b3]

[w, b]





