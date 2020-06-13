clear
load MDuchi-5class.mat;

n1 = size(X1);
n1 = n1(1);
n2 = size(X2);
n2 = n2(1);
n3 = size(X3);
n3 = n3(1);
n4 = size(X4);
n4 = n4(1);
n5 = size(X5);
n5 = n5(1);

cvx_begin
    variables w1(p) w2(p) w3(p) w4(p) w5(p) b1 b2 b3 b4 b5;
    variables slack1(n1,m-1) slack2(n2,m-1) slack3(n3,m-1) slack4(n4,m-1) slack5(n5,m-1);
    variables xi1(n1) xi2(n2) xi3(n3) xi4(n4) xi5(n5); 

    minimize( sum_square([w1;w2;w3;w4;w5])/2 + C* sum([xi1;xi2;xi3;xi4;xi5]))
    subject to
        w1 + w2 + w3 + w4 + w5 == 0;
        b1 + b2 + b3 + b4 + b5 == 0;
        X1*(w1 - w2) + b1 - b2 >= 1-slack1(:,1);
        X1*(w1 - w3) + b1 - b3 >= 1-slack1(:,2);
        X1*(w1 - w4) + b1 - b4 >= 1-slack1(:,3);
        X1*(w1 - w5) + b1 - b5 >= 1-slack1(:,4);
        X2*(w2 - w1) + b2 - b1 >= 1-slack2(:,1);
        X2*(w2 - w3) + b2 - b3 >= 1-slack2(:,2);
        X2*(w2 - w4) + b2 - b4 >= 1-slack2(:,3);
        X2*(w2 - w5) + b2 - b5 >= 1-slack2(:,4);
        X3*(w3 - w1) + b3 - b1 >= 1-slack3(:,1);
        X3*(w3 - w2) + b3 - b2 >= 1-slack3(:,2);
        X3*(w3 - w4) + b3 - b4 >= 1-slack3(:,3);
        X3*(w3 - w5) + b3 - b5 >= 1-slack3(:,4);
        X4*(w4 - w1) + b4 - b1 >= 1-slack4(:,1);
        X4*(w4 - w2) + b4 - b2 >= 1-slack4(:,2);
        X4*(w4 - w3) + b4 - b3 >= 1-slack4(:,3);
        X4*(w4 - w5) + b4 - b5 >= 1-slack4(:,4);
        X5*(w5 - w1) + b5 - b1 >= 1-slack5(:,1);
        X5*(w5 - w2) + b5 - b2 >= 1-slack5(:,2);
        X5*(w5 - w3) + b5 - b3 >= 1-slack5(:,3);
        X5*(w5 - w4) + b5 - b4 >= 1-slack5(:,4);
        slack1 >= 0;
        slack2 >= 0;
        slack3 >= 0;
        slack4 >= 0;        
        slack5 >= 0;
        xi1 >= norms_largest(slack1,1,2)/2; 
        xi2 >= norms_largest(slack2,1,2)/2;
        xi3 >= norms_largest(slack3,1,2)/2;
        xi4 >= norms_largest(slack4,1,2)/2;
        xi5 >= norms_largest(slack5,1,2)/2;
        
        xi1 >= norms_largest(slack1,2,2)/3;
        xi2 >= norms_largest(slack2,2,2)/3;
        xi3 >= norms_largest(slack3,2,2)/3;
        xi4 >= norms_largest(slack4,2,2)/3;
        xi5 >= norms_largest(slack5,2,2)/3;
        
        xi1 >= norms_largest(slack1,3,2)/4;
        xi2 >= norms_largest(slack2,3,2)/4;
        xi3 >= norms_largest(slack3,3,2)/4;
        xi4 >= norms_largest(slack4,3,2)/4;
        xi5 >= norms_largest(slack5,3,2)/4;
        
        xi1 >= norms_largest(slack1,4,2)/5;
        xi2 >= norms_largest(slack2,4,2)/5;
        xi3 >= norms_largest(slack3,4,2)/5;
        xi4 >= norms_largest(slack4,4,2)/5;
        xi5 >= norms_largest(slack5,4,2)/5;
cvx_end



[w1.' b1; w2.' b2; w3.' b3; w4.' b4; w5.' b5]




