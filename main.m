clear all;
close all;
clc;
%% System Properties and Initial codition

%Enable or Diable the Animation
animate = true;

%Central Force Constant (analogous to gravitational constant G)
G = 6.674*10^-11;

%Radius for Particles
RE = 6.4*10^6;

%Additional Particles
p1 = [0, 0, 0];
p2 = [-10 * RE + 3 * RE * rand(), 0, 0];
p3 = [-10 * RE + 3 * RE * rand(), 0, 0];
p4 = [-10 * RE + 3 * RE * rand(), 0, 0];

%Assignment of Mass to Particles
m1 = 6*10^24;
m2 = 1.5*10^4;
m3 = 1.5*10^4;
m4 = 1.5*10^4;

%Assignment of Initial Momentum to Particles
mv1 = [0 0 0];
mv2 = [0, 3*10^3 * m2, 0];
mv3 = [0, 3*10^3 * m3, 0];
mv4 = [0, 3*10^3 * m4, 0];

%Time Tnterval and Step
tmax = 100000;
dt = 100;
t = 0;

%Center of Mass
CoM = (m1*p1 + m2*p2 + m3*p3 + m4*p4)./(m1 + m2 + m3 + m4);

%Plotting the initial conditions of the simulation
figure;
hold on;
plot3(p1(1,1), p1(1,2), p1(1,3), 'pentagram');
plot3(p2(1,1), p2(1,2), p2(1,3), 'pentagram');
plot3(p3(1,1), p3(1,2), p3(1,3), 'pentagram');
plot3(p4(1,1), p4(1,2), p4(1,3), 'pentagram');
plot3(CoM(1,1), CoM(1,2), CoM(1,3), '*black');
title('Initial Conditions of the Plot');
view([45,45]);
axis equal
hold off;
%% 
% 
%% Momentum Update

%Declaration of momentum and velocity of particles


m = [m1; m2; m3; m4];
p = [p1; p2; p3; p4];



for i = 2:floor((tmax/dt)) %Run the simulation
    
    %Finding net force on each particle due to gravity
    F1 = NetForce(m, p, G, 1);
    F2 = NetForce(m, p, G, 2);
    F3 = NetForce(m, p, G, 3);
    F4 = NetForce(m, p, G, 4);
    
    %Update momentum of each particle
    mv1 = mv1 + F1 * dt;
    mv2 = mv2 + F2 * dt;
    mv3 = mv3 + F3 * dt;
    mv4 = mv4 + F4 * dt;
   
    %Computing updated position of each particle after time (dt)
    p1 = p1 + (mv1/m1) * dt;
    p2 = p2 + (mv2/m2) * dt;
    p3 = p3 + (mv3/m3) * dt;
    p4 = p4 + (mv4/m4) * dt;

    %Storing the trajectory of each particle
    p1_out(i,:) = p1;
    p2_out(i,:) = p2;
    p3_out(i,:) = p3;
    p4_out(i,:) = p4;
    p = [p1; p2; p3; p4];

    %Storing momentum as a function of time for each particle
    mv1_out(i,:) = mv1;
    mv2_out(i,:) = mv2;
    mv3_out(i,:) = mv3;
    mv4_out(i,:) = mv4;
    mv = [mv1; mv2; mv3; mv4];
    
    %Update center of mass based on new positions
    
    CoM = (m1*p1 + m2*p2 + m3*p3 + m4*p4)./(m1 + m2 + m3 + m4);

    %Storing center of mass as a function of time
    CoM_out(i,:) = CoM;
    
    %Increment the time step and log it 
    t = t + dt;
    t_out(i,:) = t;
    
end

%Plotting the trajectories of all 4 particles and the CoM

figure;
hold on;
plot3(p1_out(:,1), p1_out(:,2), p1_out(:,3), '.');
plot3(p2_out(:,1), p2_out(:,2), p2_out(:,3), '.');
plot3(p3_out(:,1), p3_out(:,2), p3_out(:,3), '.');
plot3(p4_out(:,1), p4_out(:,2), p4_out(:,3), '.');
plot3(CoM_out(2:end,1), CoM_out(2:end,2), CoM_out(2:end,3), '.black','MarkerSize',20);
view([45,45]);
xlim([0 10])
ylim([0 10])
zlim([0 10])
title('Initial Positions of particles');
hold off;


%% Animation

m = [m1; m2; m3; m4];
p = [p1; p2; p3; p4];

if animate == true
    figure;
    hold on;
    view([45,45]);
    
    for i = 1:10:length(p1_out)
        
        
        plot3(p1_out(i,1), p1_out(i,2), p1_out(i,3), '.red');
        plot3(p2_out(i,1), p2_out(i,2), p2_out(i,3), '.blue');
        plot3(p3_out(i,1), p3_out(i,2), p3_out(i,3), '.green');
        plot3(p4_out(i,1), p4_out(i,2), p4_out(i,3), '.yellow');
        plot3(CoM_out(i+1,1), CoM_out(i+1,2), CoM_out(i+1,3), '.black','MarkerSize',20);
        drawnow
        
    end
    hold off;

xlim([-57463927 100785725])
ylim([-50000000 100000000])
zlim([-1.00 1.00])
title("Trajectory Plot")
end


%% Net Force Function

function force = NetForce(m,p,G,ind)
%Let m be all masses in order m1,m2,m3,m4 (4x1)
%Let p be all positions p1,p2,p3,p4 (4x3)
%Let G be the gravitational constant
%Let ind be the index of the particle we want the net force for (values 1,2,3,4)

%Retrieve the position and mass of the indicated particle
m_ind = m(ind,1);
p_ind = p(ind,:);

%Replace the mass and position of the indicated particle with a NaN (To avoid the particle acting on itself)
m(ind, 1) = NaN;
p(ind, :) = [NaN, NaN, NaN];

%Compute the force vectors from each of the particles acting on one another

F = NaN(3,3);
count = 1;
for i = 1:length(m)
    if ~isnan(m(i,1))
        r = p(i,:) - p_ind;
        Fmag = (G*m(i,1)*m_ind)/(norm(r)^2);
        F_vec = Fmag*(r/norm(r));

        F(count,:) = F_vec;
        count = count + 1;
    end
end

%Sum the force vectors from each particle for the total net force acting on
%the indicated particle
force = sum(F);
end