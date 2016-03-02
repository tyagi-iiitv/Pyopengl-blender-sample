%Given the value of s determine the t where the value of s is achieved.
%Newton's Method is applied on B spline Cubic Curve to find the t where F(t) = g(t)-s is zero for given value s;
% s = g(t) = integration(t_min, t)(norm(y(t)))
% t(i+1) = t(i) - [F(i)/F'(i)]

clear all
close all

t_min = 0;
t_max = 1; 
dy_t  = 0;
i = 1;
tp = 0;
s = input("Enter the value of s for which you want to know\n");


%first segment of the Cublic Spline Curve
function y1 = f1 (t)
p0 = [0 0 0];
p1 = [0 0 1];
p2 = [0 1 0];
p3 = [1 0 0];
y1 = norm(((-3*t.^2 + 6*t - 3) * p0)/6 + ((9*t.^2 - 12*t) * p1)/6 + ((-9*t.^2 + 6*t + 3) * p2)/6 + (t^2* p3)/2);
endfunction

%Second segment of the Cublic Spline Curve
function y2 = f2 (t)
p1 = [0 0 1];
p2 = [0 1 0];
p3 = [1 0 0];
p4 = [1 1 1];
y2 = norm(((-3.*t.^2 + 6*t - 3) * p1)/6 + ((9.*t.^2 - 12*t) * p2)/6 + ((-9.*t.^2 + 6*t + 3) * p3)/6 + (t.^2* p4)/2);
endfunction

%Third segment of the Cublic Spline Curve
function y3 = f3 (t)
p2 = [0 1 0];
p3 = [1 0 0];
p4 = [1 1 1];
p5 = [0 1 1];
y3 = norm(((-3.*t.^2 + 6*t - 3) * p2)/6 + ((9.*t.^2 - 12*t) * p3)/6 + ((-9.*t.^2 + 6*t + 3) * p4)/6 + (t.^2* p5)/2);
endfunction

L = quad("f2",0, 1) + quad("f1",0, 1) + quad("f3",0, 1)
%intializing the value of t0
t0 = t_min + s*(t_max - t_min)/L;

if(s<L) 
%Applying Newton's Method upto i is less than 50
while (i!=50)
	
	p0 = [0 0 0];
	p1 = [0 0 1];
	p2 = [0 1 0];
	p3 = [1 0 0];
	p4 = [1 1 1];
	p5 = [0 1 1];

	%intialize the value of the time.
	t(1) = t0;

	g_t(i) = quad("f2",0, t(i)) + quad("f1",0, t(i)) + quad("f3",0, t(i));
	F = g_t(i) - s;
	dF1 = ((-3.*t(i).^2 + 6*t(i) - 3) * p0)/6 + ((9.*t(i).^2 - 12*t(i)) * p1)/6 + ((-9.*t(i).^2 + 6*t(i) + 3) * p2)/6 + (t(i).^2* p3)/2;
	dF2 = ((-3.*t(i).^2 + 6*t(i) - 3) * p1)/6 + ((9.*t(i).^2 - 12*t(i)) * p2)/6 + ((-9.*t(i).^2 + 6*t(i) + 3) * p3)/6 + (t(i).^2* p4)/2;
	dF3 = ((-3.*t(i).^2 + 6*t(i) - 3) * p2)/6 + ((9.*t(i).^2 - 12*t(i)) * p3)/6 + ((-9.*t(i).^2 + 6*t(i) + 3) * p4)/6 + (t(i).^2* p5)/2;
	tp = t(i) - (F/(norm(dF1) + norm(dF2) + norm(dF3)));
	t(i+1) = tp;
	i = i+1;
endwhile

printf("The time for given length\n%f",tp);
printf("\n");
else
 	printf("The curve length is exceeding from it's actual length.\n");
end
