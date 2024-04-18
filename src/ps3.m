%% Problem 1
IPrincipal = [7707.07451493673 0 0; ...
    0 7707.0745149367 0; ...
    0 0 18050.0227594212];
Ix = IPrincipal(1,1);
Iy = IPrincipal(2,2);
Iz = IPrincipal(3,3);

w0Deg = [8;4;6];
w0 = deg2rad(w0Deg);
tspan = 0:120;
w = eulerPropagator(w0,Ix,Iy,Iz,tspan,'Images/ps3_problem1.png');

%% Problem 2
