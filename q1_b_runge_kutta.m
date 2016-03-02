% Given the value of t determine the value of u
% 4th order Range Kutta Method 
% k1 = h∗ Sigma(t) / LengthDY ( u ) ;
% k2 = h∗ Sigma(t+h/2)/ LengthDY(u + k1 / 2) ;
% k3 = h∗ Sigma(t+h/2)/ LengthDY (u + k2 / 2) ;
% k4 = h∗ Sigma(t+h)/ LengthDY (u + k3) ;
% u =  u + k1 + 2∗( k2 + k3 ) + k4)/6 ;

clear all
close all

t_min = 0;
t_max = 1;
u_min = 0;
u_max = 1; 
dy_t  = 0;
i = 1;
tp = 0;
t = input("Enter the value of t for which you want to know u\n");

%time should be between 0 and 1
if (t<0 || t>1)
printf("The value of t should be between 0 and 1\n");
break;
endif

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
h = (t - t_min)/5000;
u0 = u_min;
t = t_min ;

%User input velocity function
function y5 = f5 (t)
y5= 2.4910;
endfunction


%Applying 4th-order Runge–Kutta method
while (i!=5000)
	 
	%intializing the initial value of u.
	u(1) = u0;
	p0 = [0 0 0];
	p1 = [0 0 1];
	p2 = [0 1 0];
	p3 = [1 0 0];
	p4 = [1 1 1];
	p5 = [0 1 1];

	% Calculating LengthDY ( u ) 
	dF1 = ((-3.*u(i).^2 + 6*u(i) - 3) * p0)/6 + ((9.*u(i).^2 - 12*u(i)) * p1)/6 + ((-9.*u(i).^2 + 6*u(i) + 3) * p2)/6 + (u(i).^2* p3)/2;
	dF2 = ((-3.*u(i).^2 + 6*u(i) - 3) * p1)/6 + ((9.*u(i).^2 - 12*u(i)) * p2)/6 + ((-9.*u(i).^2 + 6*u(i) + 3) * p3)/6 + (u(i).^2* p4)/2;
	dF3 = ((-3.*u(i).^2 + 6*u(i) - 3) * p2)/6 + ((9.*u(i).^2 - 12*u(i)) * p3)/6 + ((-9.*u(i).^2 + 6*u(i) + 3) * p4)/6 + (u(i).^2* p5)/2;

	% k1 = h∗ Sigma(t) / LengthDY ( u ) ;
	k1 = (h.*f5(t))/(norm(dF1) + norm(dF2) + norm(dF3));
	
	% Calculating the value of LengthDY(u + k1 / 2)
	dF4 = ((-3.*(u(i)+k1/2).^2 + 6*(u(i)+k1/2) - 3) * p0)/6 + ((9.*(u(i)+k1/2).^2 - 12*(u(i)+k1/2)) * p1)/6 + ((-9.*(u(i)+k1/2).^2 + 6*
(u(i)+k1/2) + 3) * p2)/6 + ((u(i)+k1/2).^2* p3)/2;
	dF5 = ((-3.*(u(i)+k1/2).^2 + 6*(u(i)+k1/2) - 3) * p1)/6 + ((9.*(u(i)+k1/2).^2 - 12*(u(i)+k1/2)) * p2)/6 + ((-9.*(u(i)+k1/2).^2 + 6*
(u(i)+k1/2) + 3) * p3)/6 + ((u(i)+k1/2).^2* p4)/2;
	dF6 = ((-3.*(u(i)+k1/2).^2 + 6*(u(i)+k1/2) - 3) * p2)/6 + ((9.*(u(i)+k1/2).^2 - 12*(u(i)+k1/2)) * p3)/6 + ((-9.*(u(i)+k1/2).^2 + 6*
(u(i)+k1/2) + 3) * p4)/6 + ((u(i)+k1/2).^2* p5)/2;

	% k2 = h∗ Sigma(t+h/2)/ LengthDY(u + k1 / 2) ;
	k2 = (h.*f5(t + h/2))/(norm(dF4) + norm(dF5) + norm(dF6));	

	% Calculating the value of LengthDY(u + k1 / 2)
	dF7 = ((-3.*(u(i)+k2/2).^2 + 6*(u(i)+k2/2) - 3) * p0)/6 + ((9.*(u(i)+k2/2).^2 - 12*(u(i)+k2/2)) * p1)/6 + ((-9.*(u(i)+k2/2).^2 + 6*
(u(i)+k2/2) + 3) * p2)/6 + ((u(i)+k2/2).^2* p3)/2;
	dF8 = ((-3.*(u(i)+k2/2).^2 + 6*(u(i)+k2/2) - 3) * p1)/6 + ((9.*(u(i)+k2/2).^2 - 12*(u(i)+k2/2)) * p2)/6 + ((-9.*(u(i)+k2/2).^2 + 6*
(u(i)+k2/2) + 3) * p3)/6 + ((u(i)+k2/2).^2* p4)/2;
	dF9 = ((-3.*(u(i)+k2/2).^2 + 6*(u(i)+k2/2) - 3) * p2)/6 + ((9.*(u(i)+k2/2).^2 - 12*(u(i)+k2/2)) * p3)/6 + ((-9.*(u(i)+k2/2).^2 + 6*
(u(i)+k2/2) + 3) * p4)/6 + ((u(i)+k2/2).^2* p5)/2;

	% k3 = h∗ Sigma(t+h/2)/ LengthDY (u + k2 / 2) ;
	k3 = (h.*f5(t + h/2))/(norm(dF7) + norm(dF8) + norm(dF9));

	% Calculating the value of LengthDY(u + k3)
	dF10 = ((-3.*(u(i)+k3).^2 + 6*(u(i)+k3) - 3) * p0)/6 + ((9.*(u(i)+k3).^2 - 12*(u(i)+k3)) * p1)/6 + ((-9.*(u(i)+k3).^2 + 6*
(u(i)+k3) + 3) * p2)/6 + ((u(i)+k3).^2* p3)/2;
	dF11 = ((-3.*(u(i)+k3).^2 + 6*(u(i)+k3) - 3) * p1)/6 + ((9.*(u(i)+k3).^2 - 12*(u(i)+k3)) * p2)/6 + ((-9.*(u(i)+k3).^2 + 6*
(u(i)+k3) + 3) * p3)/6 + ((u(i)+k3).^2* p4)/2;
	dF12 = ((-3.*(u(i)+k3).^2 + 6*(u(i)+k3) - 3) * p2)/6 + ((9.*(u(i)+k3).^2 - 12*(u(i)+k3)) * p3)/6 + ((-9.*(u(i)+k3).^2 + 6*
(u(i)+k3) + 3) * p4)/6 + ((u(i)+k3).^2* p5)/2;
	
	% k4 = h∗ Sigma(t+h)/ LengthDY (u + k3) ;
	k4  = (h.*f5(t + h))/(norm(dF10) + norm(dF11) + norm(dF12));
	
	t   = t+h; %updating the value of t
	tp = u(i) + (k1+(2*(k2+k3))+k4)/6; %updating the new value of u.
	u(i+1) = tp;
	i = i+1;
endwhile

	printf("The u for given t is %f", tp);
	printf("\n");

