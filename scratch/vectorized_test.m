function re=rabdab()
x=linspace(0,20000,20000)';
opts=odeset('Vectorized','on');

tic;
[T,Y] = ode45(@fun,[x],[0 1 1]);
[T,Y] = ode45(@fun,[x],[0 1 1]);
[T,Y] = ode45(@fun,[x],[0 1 1]);
toc;

tic;
[A,B] = ode45(@fun2,[x],[0 1 1]);
[A,B] = ode45(@fun2,[x],[0 1 1]);
[A,B] = ode45(@fun2,[x],[0 1 1]);
toc;

tic;
[A,B] = ode45(@fun3,[x],[0 1 1],opts);
[A,B] = ode45(@fun3,[x],[0 1 1],opts);
[A,B] = ode45(@fun3,[x],[0 1 1],opts);
toc;

tic;
[A,B] = ode45(@fun4,[x],[0 1 1],opts);
[A,B] = ode45(@fun4,[x],[0 1 1],opts);
[A,B] = ode45(@fun4,[x],[0 1 1],opts);
toc;


function dy = fun(t,y)
dy = zeros(3,1);    % a column vector
dy = [y(2) * y(3);...
-y(1) * y(3);...
-0.51 * y(1) * y(2);];

function dy = fun2(t,y)
dy = zeros(3,1);    % a column vector
dy(1) = y(2) * y(3);
dy(2) = -y(1) * y(3);
dy(3) = -0.51 * y(1) * y(2);

function dy = fun3(t,y)
dy = zeros(size(y)); % a matrix with arbitrarily many columns, rather than a column vector
dy = [y(2,:) .* y(3,:);...
-y(1,:) .* y(3,:);...
-0.51 .* y(1,:) .* y(2,:);];

function dy = fun4(t,y)
dy = [y(2,:) .* y(3,:);... % same as fun3()
-y(1,:) .* y(3,:);...
-0.51 .* y(1,:) .* y(2,:);];