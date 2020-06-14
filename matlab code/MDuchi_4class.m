clear
load ../MDuchi-4class.mat;

n1 = size(X1);
n1 = n1(1);
n2 = size(X2);
n2 = n2(1);
n3 = size(X3);
n3 = n3(1);
n4 = size(X4);
n4 = n4(1);
cvx_begin quiet
    variables w1(p) w2(p) w3(p) w4(p) b1 b2 b3 b4;
    variables slack1(n1,m-1) slack2(n2,m-1) slack3(n3,m-1) slack4(n4,m-1);
    variables xi1(n1) xi2(n2) xi3(n3) xi4(n4); 

    minimize( sum_square([w1;w2;w3;w4])/2 + C* sum([xi1;xi2;xi3;xi4]))
    subject to
        w1 + w2 + w3 + w4 == 0;
        b1 + b2 + b3 + b4 == 0;
        X1*(w1 - w2) + b1 - b2 >= 1-slack1(:,1);
        X1*(w1 - w3) + b1 - b3 >= 1-slack1(:,2);
        X1*(w1 - w4) + b1 - b4 >= 1-slack1(:,3);
        X2*(w2 - w1) + b2 - b1 >= 1-slack2(:,1);
        X2*(w2 - w3) + b2 - b3 >= 1-slack2(:,2);
        X2*(w2 - w4) + b2 - b4 >= 1-slack2(:,3);
        X3*(w3 - w1) + b3 - b1 >= 1-slack3(:,1);
        X3*(w3 - w2) + b3 - b2 >= 1-slack3(:,2);
        X3*(w3 - w4) + b3 - b4 >= 1-slack3(:,3);
        X4*(w4 - w1) + b4 - b1 >= 1-slack4(:,1);
        X4*(w4 - w2) + b4 - b2 >= 1-slack4(:,2);
        X4*(w4 - w3) + b4 - b3 >= 1-slack4(:,3);
        slack1 >= 0;
        slack2 >= 0;
        slack3 >= 0;
        slack4 >= 0;
        xi1 >= norms_largest(slack1,1,2)/2; 
        xi2 >= norms_largest(slack2,1,2)/2;
        xi3 >= norms_largest(slack3,1,2)/2;
        xi4 >= norms_largest(slack4,1,2)/2;
        
        xi1 >= norms_largest(slack1,2,2)/3;
        xi2 >= norms_largest(slack2,2,2)/3;
        xi3 >= norms_largest(slack3,2,2)/3;
        xi4 >= norms_largest(slack4,2,2)/3;
        
        xi1 >= norms_largest(slack1,3,2)/4;
        xi2 >= norms_largest(slack2,3,2)/4;
        xi3 >= norms_largest(slack3,3,2)/4;
        xi4 >= norms_largest(slack4,3,2)/4;
cvx_end


[w1.' b1; w2.' b2; w3.' b3; w4.' b4]


%%%%%%%%%%%%%%%
n = size(X)
n = n(1)

cvx_begin
     variables w(m,p) b(m) slack(n,m) xi(n);
    minimize( sum(sum_square(w))/2 + C* sum(xi))
    subject to
        sum(w) == 0;
        sum(b) == 0;     
        ([X, ones(n,1)]*[w,b]'.*Y_mat)*ones(m,m) - [X, ones(n,1)]*[w,b]' >= 1-slack;
        slack >= 0;
        xi >= norms_largest(slack.*(1-Y_mat),1,2)/2; 
        xi >= norms_largest(slack.*(1-Y_mat),2,2)/3
        xi >= norms_largest(slack.*(1-Y_mat),3,2)/4;
cvx_end

[w,b]

[w1.' b1; w2.' b2; w3.' b3; w4.' b4]



n = size(X)
n = n(1)

cvx_begin
     variables w(m,p) b(m) slack(n,m) xi(n);
    minimize( sum(sum_square(w))/2 + C* sum(xi))
    subject to
        sum(w) == 0;
        sum(b) == 0;     
        ([X, ones(n,1)]*[w,b]'.*Y_mat)*ones(m,m) - [X, ones(n,1)]*[w,b]' >= 1-slack;
        slack >= 0;
        xi >= norms_largest(slack.*(1-Y_mat),1,2)/2; 
        xi >= norms_largest(slack.*(1-Y_mat),2,2)/3
        xi >= norms_largest(slack.*(1-Y_mat),3,2)/4;
cvx_end

[w,b]

[w1.' b1; w2.' b2; w3.' b3; w4.' b4]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

size(permute(repmat(Y_mat,[1,1,m-1]),[1,3,2]))


cvx_begin
    variables w(m,p) b(m) slack(n,m) xi(n) t(n,m-1) u(n,m-1,m);
    minimize( sum(sum_square(w))/2 + C* sum(xi))
    subject to
        sum(w) == 0;
        sum(b) == 0;     
        ([X, ones(n,1)]*[w,b]'.*Y_mat)*ones(m,m) - [X, ones(n,1)]*[w,b]' >= 1-slack;
        slack >= 0;
        xi >= max((ones(n,1)*(1:(m-1))./(2:m)).*t  - (ones(n,1)*1./(2:m)).*sum(u,3),[],2);
        repmat(t,[1,1,m]) + u >= permute(repmat(slack,[1,1,m-1]),[1,3,2]);
        u >= 0;
%         xi >= norms_largest(slack,1,2)-1; 
%         xi >= norms_largest(slack,2,2)/2-1/2;
%         xi >= norms_largest(slack,3,2)/3-1/3;
%         xi >= 0;
cvx_end

(1:(m-1))./(2:m)
