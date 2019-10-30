%% Elastic collision between N 3D balls 

% Video: https://www.youtube.com/watch?v=78glA9qb-Pc

% DONE - There are N objects with diameters of 1 meter in a 3D space.
%        (Please define what 1 meter is in the program.)

% DONE - Each object has the same weight, but gravity can be ignored.

% DONE - Each object will start at random locations within fixed boundaries that enclose a region
%        in 3D space (a box).

% DONE - Each object will start in a random direction.

% DONE - Each object can move continously in a straight line until it comes in contact with another object or
%        it reaches a fixed boundary.

% DONE - If an object comes in contact with another object or it reaches a fixed boundary,
%   it will be reflected and change direction. Because of this, objects should not stay on the boundary
%   for a long amount of time. The objects should be reflected in a direction that is physically reasonable.

% DONE - All the objects will stay within the fixed boundaries.

% DONE - The program will terminate after a fixed amount of time.

% DONE - Save the trajectory of each object in a file.

% DONE - All the objects will move at once in one time increment.

% DONE - Write comments to explain the program.

% DONE - The movement of each object will be shown visually on the screen.

% DONE - Save snapshots of the location of the objects as image files.

% DONE - The program will terminate when the user presses Ctrl-C, but the file ...
%        with the trajectories will still be saved.

%% initials

% cleanup
clc;
clearvars;

cla;
clf;

% dock figures in predefined layout
f = figure(1);
set(f, 'DefaultFigureWindowStyle', 'docked')

fprintf('We assume that objects have same mass and perfectly elastic collision between objects and box wall,\n')
fprintf('the units of the figure are in meters and each sphere object has diameter(1m) and mass(1kg), gravity ignored\n')

% deleting old file
[st,~] = rmdir('snapshots');
if st
    mkdir snapshots
end
curDate = datestr(now,'mm-dd-yyyy HH-MM');

% check runtime
tic

% boundary limits
N=25;
M=2*N;

% 8 vertices of the boundary box
A = [N N N];
B = [N M N];
C = [N M M];
D = [N N M];
E = [M N M];
F = [M M M];
G = [M M N];
H = [M N N];

% take user input
prompt = '\nHow many objects do you wish to see (type positive integer)?\n';
objN = input(prompt);

% object properties
mass = 1;
radius = 1;
colorObj = rand(objN,3);

% object shape as sphere
t = 2*pi*rand(2,500);
P = hSphere(t,radius);
Xco=P(1,:);
Yco=P(2,:);
Zco=P(3,:);

% objects initial random position and velocities
posObj = randi([N, M],objN,3);
velObj = rand(objN,3);

%% MAIN

for i=1:200
    
    % draw plane(s) boundaries
    drawPlane(A, B, C, D, N); drawPlane(E, F, G, H, N); drawPlane(A, D, E, H, N);
    drawPlane(B, C, F, G, N); drawPlane(A, B, G, H, N); drawPlane(C, D, E, F, N);
    
    for o=1:objN
        
        % move objects in straight line with constant velocity
        posObj(o,:) = posObj(o,:) + velObj(o,:);
        
        % check boundary collision and reflect object
        for c=1:3
            
            if posObj(o,c)<=radius
                posObj(o,c) = radius;
                velObj(o,c) = -velObj(o,c);
            end
            
            if posObj(o,c)>=M-radius
                posObj(o,c) = M-radius;
                velObj(o,c) = -velObj(o,c);
            end
            
            if posObj(o,c)<=N+radius
                posObj(o,c) = N+radius;
                velObj(o,c) = -velObj(o,c);
            end
            
        end
        
        % Check for collision between each objects and reflect
        if o<objN
            
            for j=o+1:objN
                
                n = norm(posObj(o,:) - posObj(j,:));
                
                if n <= 2*radius
                    
                    % unit direction vector between object's center
                    Pnorm = (posObj(o,:) - posObj(j,:))/n;
                    
                    % initial velocities
                    u1 = velObj(o,:);
                    u2 = velObj(j,:);
                    
                    % relative velocity
                    Vrel = u1 - u2;
                    
                    % relative velocity along normal direction
                    Vnorm = dot(Vrel,Pnorm)*Pnorm;
                    
                    % since each object has same mass therefore, momentum and velocity will interchange
                    velObj(o,:) = velObj(o,:) - Vnorm;
                    velObj(j,:) = velObj(j,:) + Vnorm;
                    
                end
                
            end
            
        end
        
        % visualize
        plot3(posObj(o,1)+radius*Xco, posObj(o,2)+radius*Yco, posObj(o,3)+radius*Zco,'color',colorObj(o,:));
        title('Elastic collision between 3D objects');
        
    end
    
    hold off;
    pause(0.05); % time-step
    
    % save object Trajectories as in .mat file
    POS{i} = posObj;
    save( fullfile(pwd, strcat(curDate,'_Q1Traj.mat')), 'POS' );
    
    % save snapshots of the objects locations
    if mod(i,25)==0
        img = getframe(gcf);
        imwrite(img.cdata, [strcat('snapshots/frame_',num2str(i)), '.png'])
    end
    
end
toc

%% helper functions

function drawPlane(Pa, Pb, Pc, Pd, N)

M = 2*N;
% create box in XY plane
X = [N;M;M;N;N];
Y = [N;N;M;M;N];
Z = [0;0;0;0;0];

% create bondaries of plane
plot3(X,Y,Z+M,'r','LineWidth', 2);
plot3(X,Y,Z+N,'r','LineWidth', 2);
for j=1:length(X)-1
    plot3( [X(j);X(j)], [Y(j);Y(j)], [N;M] ,'r', 'LineWidth', 2);
end

% corners of the plane
Px = [Pa(1) Pb(1) Pc(1) Pd(1)];
Py = [Pa(2) Pb(2) Pc(2) Pd(2)];
Pz = [Pa(3) Pb(3) Pc(3) Pd(3)];

patch(Px, Py, Pz, 'y');
alpha(0.1);

% maintain axis limits
scal = 3;
xlim([N-scal M+scal]); xlabel('X(m)');
ylim([N-scal M+scal]); ylabel('Y(m)');
zlim([N-scal M+scal]); zlabel('Z(m)');
axis square manual; grid on; box on;

hold on;

end

function S = hSphere(sample, r)

rows = size(sample,1);
cols = size(sample,2);
out = ones(rows+1, cols);

for i = 1:rows
    j = repmat(cos(sample(i, :)), i, 1);
    j = [j; sin(sample(i, :))];
    j = [j; ones(rows - i, cols)];
    out = out.*j;
end
S = out*r;

end

