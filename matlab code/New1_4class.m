clear
load New1-4class.mat;

n1 = size(X1);
n1 = n1(1);
n2 = size(X2);
n2 = n2(1);
n3 = size(X3);
n3 = n3(1);
n4 = size(X4);
n4 = n4(1);
cvx_begin quiet
    variables w1(p) w2(p) w3(p) b1 b2 b3;
    variables slack1(n1,m-1) slack2(n2,m-1) slack3(n3,m-1) slack4(n4,m-1);
    variables xi1(n1) xi2(n2) xi3(n3) xi4(n4); 

    minimize( sum_square([w1;w2;w3])/2 + C* (sum([xi1;xi2;xi3]) + sum(sum(slack4))))
    subject to
        slack1(:,1) == 1;
        slack2(:,2) == 1;
        slack3(:,3) == 1;
        X1*(w1 - w2) + b1 - b2 >= 1-slack1(:,2);
        X1*(w1 - w3) + b1 - b3 >= 1-slack1(:,3);
        X1*w1 + b1 >= 1 - xi1;
        X2*(w2 - w1) + b2 - b1 >= 1-slack2(:,1);
        X2*(w2 - w3) + b2 - b3 >= 1-slack2(:,3);
        X2*w2 + b2 >= 1 - xi2;
        X3*(w3 - w1) + b3 - b1 >= 1-slack3(:,1);
        X3*(w3 - w2) + b3 - b2 >= 1-slack3(:,2);
        X3*w3 + b3 >= 1 - xi3;
        X4*w1 + b1 <= -0 + slack4(:,1);
        X4*w2 + b2 <= -0 + slack4(:,2);
        X4*w3 + b3 <= -0 + slack4(:,3);
        slack1 >= 0;
        slack2 >= 0;
        slack3 >= 0;
        slack4 >= 0;
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
