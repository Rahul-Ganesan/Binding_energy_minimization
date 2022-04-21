function min = CnG(f)
syms b phi theta
% initial guess vaules
x(1) = 1000;
y(1) = -200;
z(1) = 50;

% gradient computation
df_dx = diff(f, b);
df_dy = diff(f, phi);
df_dz = diff(f, theta);
g0 = [subs(df_dx, [b,phi,theta], [x(1),y(1),z(1)]) 
      subs(df_dy, [b,phi,theta], [x(1),y(1),z(1)]) 
      subs(df_dz, [b,phi,theta], [x(1),y(1),z(1)]) ];

% search direction 
d = -g0;

% iterations
for i = 1:20
    I = [x(i),y(i),z(i)];
    syms t 
    g = subs(f, [b,phi,theta], [x(i)+d(1)*t,y(i)+t*d(2),z(i)+t*d(3)]);
    dg_dt = diff(g,t);
    t = solve(dg_dt,t);
    x(i+1) = I(1)+t*d(1); % New x value
    y(i+1) = I(2)+t*d(2); % New y value
    z(i+1) = I(3)+t*d(3); % New z value
    g_o = [subs(df_dx, [b,phi,theta], [x(i),y(i),z(i)]) 
           subs(df_dy, [b,phi,theta], [x(i),y(i),z(i)])
           subs(df_dz, [b,phi,theta], [x(i),y(i),z(i)])]
    i = i+1;
    g_1 = [subs(df_dx, [b,phi,theta], [x(i),y(i),z(i)]) 
           subs(df_dy, [b,phi,theta], [x(i),y(i),z(i)])
           subs(df_dz, [b,phi,theta], [x(i),y(i),z(i)])]  % Updated Gradient
    beta = (g_1/g_o)^2
    d
    d=-(g_1)+beta*d;
end
Iter = 1:i;
X_coordinate = x';
Y_coordinate = y';
Z_coordinate = z';
Iterations = Iter';
T = table(Iterations,X_coordinate,Y_coordinate,Z_coordinate)

min= [ x(i), y(i), z(i), subs(f,[b,phi,theta], [x(i),y(i),z(i)])];
end