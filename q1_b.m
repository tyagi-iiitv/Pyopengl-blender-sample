% Given the value of s determine the t where the value of s is achieved.
% Newton's Method is applied on B spline Cubic Curve to find the u for given value of t;
% Now s is varying entity. We have speed function that is varying with t.
% s = integration(0, given value of t)(speed_function)
% s = g(t) = integration(t_min, t)(norm(y(t)))
% t(i+1) = t(i) - [F(i)/F'(i)]

clear all
close all
s = 0;
t_min = 0;
t_max = 1; 
u_min = 0;
u_max = 1;
dy_t  = 0;
i = 1;
tp = 0;
t = input("Enter the value of t for which you want to know value of u\n");


%User specified speed function
function y4 = f4 (t)
y4 = 2.4910;
endfunction

%first segment of the Cublic Spline Curve
function y1 = f1 (u)
p0 = [0 0 0];
p1 = [0 0 1];
p2 = [0 1 0];
p3 = [1 0 0];
y1 = norm(((-3*u.^2 + 6*u - 3) * p0)/6 + ((9*u.^2 - 12*u) * p1)/6 + ((-9*u.^2 + 6*u + 3) * p2)/6 + (u^2* p3)/2);
endfunction

%Second segment of the Cublic Spline Curve
function y2 = f2 (u)
p1 = [0 0 1];
p2 = [0 1 0];
p3 = [1 0 0];
p4 = [1 1 1];
y2 = norm(((-3.*u.^2 + 6*u - 3) * p1)/6 + ((9.*u.^2 - 12*u) * p2)/6 + ((-9.*u.^2 + 6*u + 3) * p3)/6 + (u.^2* p4)/2);
endfunction

%Third segment of the Cublic Spline Curve
function y3 = f3 (u)
p2 = [0 1 0];
p3 = [1 0 0];
p4 = [1 1 1];
p5 = [0 1 1];
y3 = norm(((-3.*u.^2 + 6*u - 3) * p2)/6 + ((9.*u.^2 - 12*u) * p3)/6 + ((-9.*u.^2 + 6*u + 3) * p4)/6 + (u.^2* p5)/2);
endfunction

L = quad("f2",0, u_max) + quad("f1",0, u_max) + quad("f3",0, u_max)
%intializing the value of t0

if ( t>=0 && t<=1)
s = quad("f4",0, t);
u0 = u_min + s*(u_max - u_min)/L;
%Applying Newton's Method upto i is less than 50
while (i!=50)
	
	p0 = [0 0 0];
	p1 = [0 0 1];
	p2 = [0 1 0];
	p3 = [1 0 0];
	p4 = [1 1 1];
	p5 = [0 1 1];

	%intialize the value of the time.
	u(1) = u0;

	g_t(i) = quad("f2",0, u(i)) + quad("f1",0, u(i)) + quad("f3",0, u(i));
	F = g_t(i) - s;
	dF1 = ((-3.*u(i).^2 + 6*u(i) - 3) * p0)/6 + ((9.*u(i).^2 - 12*u(i)) * p1)/6 + ((-9.*u(i).^2 + 6*u(i) + 3) * p2)/6 + (u(i).^2* p3)/2;
	dF2 = ((-3.*u(i).^2 + 6*u(i) - 3) * p1)/6 + ((9.*u(i).^2 - 12*u(i)) * p2)/6 + ((-9.*u(i).^2 + 6*u(i) + 3) * p3)/6 + (u(i).^2* p4)/2;
	dF3 = ((-3.*u(i).^2 + 6*u(i) - 3) * p2)/6 + ((9.*u(i).^2 - 12*u(i)) * p3)/6 + ((-9.*u(i).^2 + 6*u(i) + 3) * p4)/6 + (u(i).^2* p5)/2;
	tp = u(i) - (F/(norm(dF1) + norm(dF2) + norm(dF3)));
	u(i+1) = tp;
	i = i+1;
endwhile

printf("The value of u for the given value of t is\n%f",tp);
printf("\n");
else
	printf("The time entered should be between 0 and 1.\n");
end
