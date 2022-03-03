Closing IO
stagesT.moveAbsAx('X', 13)
x=0:0.1:100;
y=exp(-x/10);
figure; plot(x,y)
n=randn(1,1000)*0.05;
figure; plot(n)
yy=sqrt(exp(-x/10)^2+n.^2);
Error using  ^  (line 51)
Incorrect dimensions for raising a matrix to a power. Check that the matrix is
square and the power is a scalar. To perform elementwise matrix powers, use '.^'.
 
yy=sqrt(exp(-x/10)^2+n'.^2);
Error using  ^  (line 51)
Incorrect dimensions for raising a matrix to a power. Check that the matrix is
square and the power is a scalar. To perform elementwise matrix powers, use '.^'.
 
yy=sqrt(exp(-x/10).^2+n.^2);
Matrix dimensions must agree.
 
yy=sqrt(exp(-x'/10).^2+n.^2);
figure; plot(y)
figure; plot(yy)
size(x)

ans =

           1        1001

size(n)

ans =

           1        1000

n=randn(1,1001)*0.05;
yy=sqrt(exp(-x'/10).^2+n.^2);
yy=sqrt(exp(-x/10).^2+n.^2);
figure; plot(yy)
figure; plot(log(yy))
figure; plot(n)
yy=0;
for 1:100
 for 1:100
     ↑
Error: Invalid expression. Check for missing multiplication operator, missing or
unbalanced delimiters, or other syntax error. To construct matrices, use brackets
instead of parentheses.
 
for ii=1:100
yy=yy+sqrt(exp(-x'/10).^2+n.^2);
end
 for ii=1:100
n=randn(1,1001)*0.05;
yy=yy+sqrt(exp(-x'/10).^2+n.^2);
end
yy=0

yy =

     0

for ii=1:100
n=randn(1,1001)*0.05;
yy=yy+sqrt(exp(-x'/10).^2+n.^2);
end
figure; plot(y.^2)
figure; plot(log(y))
figure; plot(log(yy))
 yy=0

yy =

     0

for ii=1:100
n=randn(1,1001)*0.05;
yy=yy+sqrt(exp(-x/10).^2+n.^2);
end
figure; plot(log(yy))
for ii=1:1000
n=randn(1,1001)*0.05;
yy=yy+sqrt(exp(-x/10).^2+n.^2);
end
figure; plot(log(yy))
figure; plot((yy))
Y=sqrt(yy.^2-44^2);
figure; plot(Y)
Y=yy.^2-44^2;
figure; plot(yy.^2)
Y=sqrt(abs(yy.^2-1950));
figure; plot(Y)
figure; plot(log(Y))
hold on
plot(log(Y))
for ii=1:10000
n=randn(1,1001)*0.05;
yy=yy+sqrt(exp(-x/10).^2+n.^2);
end
figure; plot(yy.^2)
Y=sqrt(abs(yy.^2-19500));
figure; plot(log(Y))
figure; plot(Y)
figure; plot(yy)
clear yy Y x n
x=0:0.1:100;
y=0

y =

     0

yy=0

yy =

     0

 for ii=1:10000
n=randn(1,1001)*0.05;
yy=yy+sqrt(exp(-x/10).^2+n.^2);
end
figure; plot(yy.^2)
Y=sqrt(abs(yy.^2-1.58e5));
figure; plot(Y)
figure; plot(log(Y))
hold on
plot(log(Y/10))
