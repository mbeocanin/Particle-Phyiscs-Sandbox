% MIT License

% Copyright (c) 2022 Miloš Beočanin

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:

% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

% Example 4: Rope
% Demonstrates physics simulation model presented in this article:
% <https://box2d.org/files/ErinCatto_IterativeDynamics_GDC2005.pdf>

clc
clear all

% simulation constants
fps = 60;
timeScale = 1.0;

% constraint solver constants
bias = 0.2 % maximum penetration error corrective impulse per simulation iteration; lower values improve stability, but increase separation time
slop = 0.005 % allowed penetration; lower values increase separation, but reduce stability
iterations = 10 % maximum constraint resolution iterations per simulation iteration; lower values improve performance, but reduce stability
warmStart = true % solution caching between simulation iterations; improves stability while allowing for less constraint resolution iterations

% physics constants
worldSize = [2 2]; % [m]

g = 9.807; % [m/s^2]; gravitational acceleration
airDensity = 1.225; % [kg/m^3]
hempDensity = 860; % [kg/m^3]
dragCoefficient = 0.47; % for spheres

% set up objects
circleCount = 15;

% [edge1 edge2 edge3 edge4 fixedCircle circle2 ... circleN]
edgeCount = 4;
objectCount = edgeCount + circleCount; % 4 edges + circles
r = zeros(objectCount, 1); % [m]; radii
m = Inf(objectCount, 1); % [kg]; masses
v = zeros(objectCount, 2); % [m/s]; velocities
p = zeros(objectCount, 2); % [m]; positions

A = zeros(objectCount, 1); % reference areas
Fweight = zeros(objectCount, 2); % [N = kg*m/s^2]; weights

% set up edges
p(1, :) = [worldSize(1) worldSize(2)]; % top
p(2, :) = [0 0]; % bottom
p(3, :) = [0 worldSize(2)]; % left
p(4, :) = [worldSize(1) 0]; % right

normalEdge = zeros(4, 2); % edge normals
normalEdge(1, :) = [0 -1]; % top
normalEdge(2, :) = [0 1]; % bottom
normalEdge(3, :) = [1 0]; % left
normalEdge(4, :) = [-1 0]; % right

% set up rope vertices
circleRadius = 0.025; % [m]
r(edgeCount + 1:end) = ones(circleCount, 1)*circleRadius; % skip edges
m(edgeCount + 2:end) = hempDensity*4/3*r(edgeCount + 2:end).^3*pi; % skip edges, fixed vertex
for circle = 1:circleCount
    p(edgeCount + circle, :) = [worldSize(1)*0.5 worldSize(2) - worldSize(2)*0.75/circleCount*(circle - 1)]; % skip edges
end

A(edgeCount + 1:end) = r(edgeCount + 1:end).^2*pi; % reference areas for spheres, skip edges
Fweight(edgeCount + 2:end, :) = [zeros(circleCount - 1, 1) -m(edgeCount + 2:end)*g]; % skip edges, fixed vertex

% set up rope segments
segmentLength = 0.1; % [m]

% set up user force
rUser = min(worldSize)*0.1; % [m]
FUser = 5.0; % [N = kg*m/s^2]
pUser = [-Inf -Inf]; % [m]

% precomputed values
mInv = 1.0./m; % [1/kg]; inverse masses
MInv = diag(repelem(mInv, 2)); % diagonal matrix of inverse masses

% set up GUI
fig = figure('Name', 'Example 4: Rope', 'Units', 'normalized', 'Position', [0.2 0.2 0.6 0.6]);
constraintSolverPanel = uipanel(fig, 'Title', 'Constraint solver', 'Position', [0.1 0.1 0.8 0.14]);
uicontrol('Parent', constraintSolverPanel, 'Style', 'text', 'String', 'bias:', 'HorizontalAlignment', 'right', 'Units', 'normalized', 'Position', [0 0.76 0.19 0.20])
uicontrol('Parent', constraintSolverPanel, 'Style', 'text', 'String', 'slop:', 'HorizontalAlignment', 'right', 'Units', 'normalized', 'Position', [0 0.51 0.19 0.20])
uicontrol('Parent', constraintSolverPanel, 'Style', 'text', 'String', 'iterations:', 'HorizontalAlignment', 'right', 'Units', 'normalized', 'Position', [0 0.26 0.19 0.20])
uicontrol('Parent', constraintSolverPanel, 'Style', 'text', 'String', 'warm start:', 'HorizontalAlignment', 'right', 'Units', 'normalized', 'Position', [0 0.01 0.19 0.20])
uicontrol('Parent', constraintSolverPanel, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.2 0.76 0.8 0.24], 'Min', 0.0, 'Max', 1.0, 'SliderStep', [0.001 0.01], 'Value', bias, 'Callback', @guiBiasSlider);
uicontrol('Parent', constraintSolverPanel, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.2 0.51 0.8 0.24], 'Min', 0.0, 'Max', 0.1, 'SliderStep', [0.01 0.1], 'Value', slop, 'Callback', @guiSlopSlider);
uicontrol('Parent', constraintSolverPanel, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.2 0.26 0.8 0.24], 'Min', 1.0, 'Max', 100.0, 'SliderStep', [0.01 0.1], 'Value', iterations, 'Callback', @guiIterationsSlider);
uicontrol('Parent', constraintSolverPanel, 'Style', 'checkbox', 'Units', 'normalized', 'Position', [0.2 0.01 0.2 0.24], 'Value', warmStart, 'Callback', @guiWarmStartCheckbox); 
simulationPanel = uipanel(fig, 'Title', 'Simulation', 'Position', [0.1 0.01 0.8 0.08]);
uicontrol('Parent', simulationPanel, 'Style', 'text', 'String', 'time scale:', 'HorizontalAlignment', 'right', 'Units', 'normalized', 'Position', [0 0.51 0.19 0.39])
uicontrol('Parent', simulationPanel, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.2 0.51 0.8 0.49], 'Min', 0.1, 'Max', 10.0, 'SliderStep', [0.01 0.1], 'Value', timeScale, 'Callback', @guiTimeScaleSlider);
uicontrol('Parent', simulationPanel, 'Style', 'togglebutton', 'String', 'Pause', 'Units', 'normalized', 'Position', [0.2 0.01 0.1 0.49], 'Callback', @guiPauseTogglebutton);
uicontrol('Parent', simulationPanel, 'Style', 'pushbutton', 'String', 'Quit', 'Units', 'normalized', 'Position', [0.31 0.01 0.1 0.49], 'Callback', @guiQuitPushbutton);

set(fig, 'WindowButtonMotionFcn', @guiMouseMove);

% axes
set(gca, 'OuterPosition', [0.1 0.2 0.8 0.8])
axis([0 worldSize(1) 0 worldSize(2)])
axis square
axis off

% global variables
setappdata(fig, 'bias', bias)
setappdata(fig, 'slop', slop)
setappdata(fig, 'iterations', iterations)
setappdata(fig, 'warmStart', warmStart)

setappdata(fig, 'running', true)
setappdata(fig, 'paused', false)
setappdata(fig, 'timeScale', timeScale)

setappdata(fig, 'forceEmitterLocation', pUser)

% graphics objects
% [edge1 edge2 edge3 edge4 ropeSegment1 ropeSegment2 ... ropeSegmentN ropeVertex1 ropeVertex2 ... ropeVertexN]
gObjectCount = objectCount + circleCount - 1;
graphics = gobjects(gObjectCount);
% edges
graphics(1) = line([0 worldSize(1)], [worldSize(2) worldSize(2)], 'color', 'black'); % top
graphics(2) = line([0 worldSize(1)], [0 0], 'color', 'black'); % bottom
graphics(3) = line([0 0], [0 worldSize(2)], 'color', 'black'); % left
graphics(4) = line([worldSize(1) worldSize(1)], [0 worldSize(2)], 'color', 'black'); % right
% rope segments
for segment = objectCount + 1:gObjectCount % skip edges, skip vertices
    % first vertex position
    x1 = p(segment - circleCount, 1);
    y1 = p(segment - circleCount, 2);
    % second vertex position
    x2 = p(segment - circleCount + 1, 1);
    y2 = p(segment - circleCount + 1, 2);

    graphics(segment) = line([x1 x2], [y1 y2], 'LineWidth', 10.0, 'Color', 'red');
end
% rope vertices
for circle = edgeCount + 1:objectCount % skip edges
	radius = r(circle);
	diameter = 2*radius;
	x = p(circle, 1) - radius;
	y = p(circle, 2) - radius;

	position = [x y diameter diameter];
    graphics(circle) = rectangle('Position', position, 'Curvature', [1 1], 'EdgeColor', 'red', 'FaceColor', 'red'); % circles are drawn as squares with curved corners
end

% run simulation
oldConstraints = [];

while ishandle(fig)
    % read global variables
    bias = getappdata(fig, 'bias');
    slop = getappdata(fig, 'slop');
    iterations = getappdata(fig, 'iterations');
    warmStart = getappdata(fig, 'warmStart');

    running = getappdata(fig, 'running');
    paused = getappdata(fig, 'paused');
    timeScale = getappdata(fig, 'timeScale');

    pUser = getappdata(fig, 'forceEmitterLocation'); 

    % simulation control
    if ~running
        break
    elseif paused
        pause(1/fps);
        continue
    end

    % apply forces
    % -----------------------------------------------------------------------------------------------------------------
    Fext = zeros(objectCount, 2); % external forces
    for object = 1:objectCount
        % aerodynamic drag
        velocity = v(object, :);
        Fdrag = -velocity*norm(velocity)*0.5*airDensity*dragCoefficient*A(object);

        % user force
        Fuser = [0 0];
        direction = p(object, :) - pUser(1, 1:2);
        if norm(direction) <= rUser
            Fuser = normalize(direction)*FUser;
        end

        % total external force
        Fext(object, :) = Fweight(object, :) + Fdrag + Fuser;
    end
    
    dt = 1.0/fps*timeScale;
    
    % constraint setup
    % -----------------------------------------------------------------------------------------------------------------
    constraints = [];
    
    % collision detection; O(n^2) broad phase/narrow phase
    uniquePairs = nchoosek(1:objectCount, 2); % skip self contact, skip duplicate contacts
    for pair = 1:size(uniquePairs, 1);
        object1 = uniquePairs(pair, 1);
        object2 = uniquePairs(pair, 2);

        % skip immovable objects contact
        if mInv(object1) == 0 && mInv(object2) == 0
            continue
        end
        
        if object1 <= edgeCount
            % first body is edge
            normal = normalEdge(object1, :);

            direction = p(object2, :) - p(object1, :);
            separation = dot(direction, normal) - r(object2);
            if separation > 0
                continue
            end
        elseif object2 <= edgeCount
            % second body is edge
            normal = normalEdge(object2, :);

            direction = p(object1, :) - p(object2, :);
            separation = dot(direction, normal) - r(object1);
            if separation > 0
                continue
            end
        else
            % both bodies are circles
            direction = p(object2, :) - p(object1, :);
            normal = normalize(direction);
                
            separation = norm(direction) - r(object1) - r(object2);
            if separation > 0
                continue
            end
        end

        separation = min(separation + slop, 0); % allow a bit of penetration (slop); must remain negative value
        initialSolutionVelocity = 0;
        initialSolutionPosition = 0;
        solutionLowerLimit = 0; % for contact constraints
        solutionUpperLimit = Inf; % for contact constraints
        constraints = [constraints; object1 object2 normal(1) normal(2) separation initialSolutionVelocity initialSolutionPosition solutionLowerLimit solutionUpperLimit];
    end
    % rope segments
    for segment = objectCount + 1:objectCount + circleCount - 1 % skip edges, circles
        object1 = segment - circleCount;  % first vertex
        object2 = segment - circleCount + 1; % second vertex

        direction = p(object2, :) - p(object1, :);
        normal = normalize(direction);

        separation = norm(direction) - segmentLength;

        initialSolutionVelocity = 0;
        initialSolutionPosition = 0;
        solutionLowerLimit = -Inf; % for distance constraints
        solutionUpperLimit = Inf; % for distance constraints
        constraints = [constraints; object1 object2 normal(1) normal(2) separation initialSolutionVelocity initialSolutionPosition solutionLowerLimit solutionUpperLimit];
    end
    
	% constraint solving
    % -----------------------------------------------------------------------------------------------------------------
    FcVelocity = zeros(objectCount, 1); % velocity constraint forces
    FcPosition = zeros(objectCount, 1); % position constraint forces

    constraintCount = size(constraints, 1);
    if constraintCount > 0
        V = reshape(v', [], 1); % single column vector: [x1 y1 x2 y2 x3 y3 ...]'
        FExt = reshape(Fext', [], 1); % single column vector: [x1 y1 x2 y2 x3 y3 ...]'

        if warmStart
            % restore old lambda from cache
            for constraint = 1:constraintCount
                object1 = constraints(constraint, 1);
                object2 = constraints(constraint, 2);
                for oldContact = 1:size(oldConstraints, 1)
                    oldObject1 = oldConstraints(oldContact, 1);
                    oldObject2 = oldConstraints(oldContact, 2);
                    if (oldObject1 == object1 && oldObject2 == object2) || (oldObject1 == object2 && oldObject2 == object1)
                        constraints(constraint, 6) = oldConstraints(oldContact, 6);
                        %constraints(constraint, 7) = oldConstraints(oldContact, 7); % position magnitudes should not be cached
                        break
                    end
                end
            end
        end

        J = zeros(constraintCount, objectCount*2); % Jacobian
        S = zeros(constraintCount, 1); % constraint violation measure
        lambda0Velocity = zeros(constraintCount, 1); % velocity constraint initial solution
        lambda0Position = zeros(constraintCount, 1); % position constraint initial solution
        lambdaMin = -Inf(constraintCount, 1); % solution lower limit
        lambdaMax = Inf(constraintCount, 1); % solution upper limit
        for constraint = 1:constraintCount
            object1 = constraints(constraint, 1);
            object2 = constraints(constraint, 2);
            normal = constraints(constraint, 3:4);

            S(constraint, 1) = constraints(constraint, 5);
            lambda0Velocity(constraint, 1) = constraints(constraint, 6);
            lambda0Position(constraint, 1) = constraints(constraint, 7);
            lambdaMin(constraint, 1) = constraints(constraint, 8);
            lambdaMax(constraint, 1) = constraints(constraint, 9);

            % -object1 + object2
            J(constraint, object1*2 - 1) = -normal(1); % x
            J(constraint, object1*2) = -normal(2); % y
            J(constraint, object2*2 - 1) = normal(1); % x
            J(constraint, object2*2) = normal(2); % y
        end

        % MLCP (Mixed Linear Complementarity Problem)
        left = J*MInv*J'; % effective mass inverse
        rightVelocity = -J*(1/dt*V + MInv*FExt); % -(velocity acceleration + external force acceleration)      
        rightPosition = -bias/dt^2*S; % -constraint violation acceleration 

        lambdaVelocity = projectedGS(left, rightVelocity, lambda0Velocity, lambdaMin, lambdaMax, iterations, 10^-4); % velocity force magnitude
        FcVelocity = J'*lambdaVelocity; % velocity constraint force; single column vector: [x1 y1 x2 y2 x3 y3...]'

        lambdaPosition = projectedGS(left, rightPosition, lambda0Position, lambdaMin, lambdaMax, iterations, 10^-4); % position force magnitude
        FcPosition = J'*lambdaPosition; % position constraint force; single column vector: [x1 y1 x2 y2 x3 y3...]'
        
        if warmStart
            % cache lambda
            constraints(:, 6) = lambdaVelocity;
            %constraints(:, 7) = lambdaPosition; % position magnitudes should not be cached
            oldConstraints = constraints;
        end
     
        FcVelocity = reshape(FcVelocity, 2, [])'; % two column matrix [x1 y1; x2 y2; x3 y3; ...]
        FcPosition = reshape(FcPosition, 2, [])'; % two column matrix [x1 y1; x2 y2; x3 y3; ...]
    end

    % integration
    % -----------------------------------------------------------------------------------------------------------------
    for object = 1:objectCount
        % 2nd Newton's Law (a = F/m)
        aVelocity = (Fext(object, :) + FcVelocity(object, :))*mInv(object); % velocity acceleration
        aPosition = FcPosition(object, :)*mInv(object); % position acceleration

        % Semi-implicit Euler Integrator
        v(object, :) = v(object, :) + aVelocity*dt;
        vPosition = aPosition*dt;
        
        p(object, :) = p(object, :) + (v(object, :) + vPosition)*dt;
    end

    % rendering; update positions of already created graphics objects
    % -----------------------------------------------------------------------------------------------------------------
    for segment = objectCount + 1:gObjectCount % skip edges, skip vertices
        x1 = p(segment - circleCount, 1);
        y1 = p(segment - circleCount, 2);
        x2 = p(segment - circleCount + 1, 1);
        y2 = p(segment - circleCount + 1, 2);

        set(graphics(segment), 'x', [x1 x2]);
        set(graphics(segment), 'y', [y1 y2]);
    end
    for circle = edgeCount + 1:objectCount % skip edges
        radius = r(circle);
        diameter = 2*radius;
        x = p(circle, 1) - radius;
        y = p(circle, 2) - radius;

        position = [x y diameter diameter];
        set(graphics(circle), 'Position', position);
    end
    
    pause(1/fps);
end

close