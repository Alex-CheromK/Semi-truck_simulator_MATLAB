%% Setup workspace
clear;
clc;

%% Define simulation settings
% Simulation duration (sec)
Sim.EndTime = 300;
% Time interval for saving results (sec).
% Note: Sim.SavePeriod must be at most half of the controller sampling
% period. It will automatically be reduced later on.
Sim.SavePeriod = 0.1;
% Compute the number of saved samples
Sim.NumSavePoints = Sim.EndTime/Sim.SavePeriod + 1;
% Playback rate of the animation. Example: Enter "2" to play the animation
% at twice the rate of the simulated physics.
% Enter "0" for no animation.
Sim.AnimPlayRate = 1;

%% Define controller parameters
% Controller sampling period (sec)
Control.SamplingPeriod = 10;
% Compute number of controller sample points
Control.NumSamples = Sim.EndTime/Control.SamplingPeriod + 1;
% Initial control values
Control.ThrustForce = 0;
Control.WheelTorque = 0;
% Set Sim.SavePeriod to at most half of Control.SamplingPeriod
Sim.SavePeriod = min(Sim.SavePeriod,Control.SamplingPeriod/2);

%% Define global parameters
% Gravitational acceleration (m/s^2)
Global.GravAccel = 9.81;
% Rolling friction coefficient between tires and road
Global.DynFricCoeff = 0.002;
% Density of air (kg/m^3)
Global.AirDensity = 1.2;

%% Define truck properties
% Mass (kg)
Truck.Mass = 1500;
% Moment of inertia about vertical axis (kg-m^2)
Truck.MOI = 10000;
% Length from CG to front wheels (m)
Truck.CG2FrontLength = 1;
% Length from CG to rear wheels (m)
Truck.CG2RearLength = 1;
% Front wheel system moment of inertia about vertical axis (kg-m^2)
Truck.WheelMOI = 1;
% Gear ratio from steering wheel to front wheels
Truck.WheelGearRatio = 1;
% Wheel system linear rotational stiffness coefficient (N/rad)
Truck.WheelStiffCoeff = 0;
% Wheel system pivot offset due to caster angle (m)
Truck.WheelOffset = 0.25;
% Wheel system linear rotational damping coefficient (N/(rad/s))
Truck.WheelDampCoeff = 10;
% Truck-trailer system aerodynamic drag coefficient
Truck.DragCoeff = 5;
% Truck-trailer system reference drag area (m^2)
Truck.DragArea = 3;
% Artifical maximum torque for limiting wheel rotation beyond maximum
% allowable value (N-m)
Truck.WheelStopStiff = 100;
% Maximum allowable wheel angle (deg)
Truck.MaxWheelAngle = 10;

%% Define trailer properties
% Mass (kg)
Trailer.Mass = 3000;
% Moment of inertia about vertical axis (kg-m^2)
Trailer.MOI = 20000;
% Length from CG to pivot point (m)
Trailer.Pivot2CG_Length = 6;
% Length from CG to rear wheels (m)
Trailer.Pivot2RearLength = 12;

%% Initialize states and results class objects
% Initialize state vector
InitStateVec = [
    0; % Truck x position in global frame (m)
    0; % Truck y position in global frame (m)
    45*pi/180; % Truck angle (rad)
    0; % Truck forward velocity in truck frame (m/s)
    0; % Trailer angle relative to truck (rad)
    0; % Front wheel angle (rad)
    0]; % Front wheel angular velocity (rad/s)

% Initialize the simulation start time (sec)
StartTime = 0;

% Save simulation results at time = zero
for i = Sim.NumSavePoints:-1:1
    Results.TimeStep(i).Time = StartTime;
    Results.TimeStep(i).StateVec = InitStateVec;
    [~,...
        Results.TimeStep(i).Truck,...
        Results.TimeStep(i).Trailer] = TruckDynModel(...
        StartTime,...
        InitStateVec,...
        Control,Global,Truck,Trailer);
end
TimeStepNum = 2;

% Save results at the first controller sample
for i = Control.NumSamples:-1:1
    Results.Sample(i).Time = StartTime;
    Results.Sample(i).StateVec = InitStateVec;
    [~,...
        Results.Sample(i).Truck,...
        Results.Sample(i).Trailer] = TruckDynModel(...
        StartTime,...
        InitStateVec,...
        Control,Global,Truck,Trailer);
end

%% Run discrete-time simulation
for SampleNum = 1:Control.NumSamples - 1
    % Define control inputs for the current sampling period
    Control.ThrustForce = 100; % N
    Control.WheelTorque = 1; % N-m
    % Note: When testing controllers, the above terms may be calculated
    % from a control function.
    
    % Solve dynamics problem for the current smapling period
    [Sol.TimeTraj,Sol.StateVecTraj] = ode45(...
        @(Time,StateVec) TruckDynModel(...
        Time,StateVec,Control,Global,Truck,Trailer),...
        (StartTime:Sim.SavePeriod:StartTime + Control.SamplingPeriod),...
        InitStateVec);
    
    % Save simulation results
    for i = 2:length(Sol.TimeTraj)
        Results.TimeStep(TimeStepNum).Time = Sol.TimeTraj(i);
        Results.TimeStep(TimeStepNum).StateVec = Sol.StateVecTraj(i,:)';
        [~,...
            Results.TimeStep(TimeStepNum).Truck,...
            Results.TimeStep(TimeStepNum).Trailer] = TruckDynModel(...
            Results.TimeStep(TimeStepNum).Time,...
            Results.TimeStep(TimeStepNum).StateVec,...
            Control,Global,Truck,Trailer);
        TimeStepNum = TimeStepNum + 1;
    end
    
    % Save results at the end of the current controller sampling period
    Results.Sample(SampleNum + 1).Time = Sol.TimeTraj(end);
    Results.Sample(SampleNum + 1).StateVec = Sol.StateVecTraj(end,:)';
    [~,...
        Results.Sample(SampleNum + 1).Truck,...
        Results.Sample(SampleNum + 1).Trailer] = TruckDynModel(...
        Results.Sample(SampleNum + 1).Time,...
        Results.Sample(SampleNum + 1).StateVec,...
        Control,Global,Truck,Trailer);
    
    % Update ODE solver inputs for the next sampling period
    StartTime = StartTime + Control.SamplingPeriod;
    InitStateVec = Sol.StateVecTraj(end,:)';
end

%% Generate truck-trailer animation
if Sim.AnimPlayRate > 0
    % Define truck dimensions for visualization purposes only (m)
    Truck.Width = 2;
    Trailer.Width = 2.5;
    
    % Initialize a points vector for plotting purposes
    p = zeros(2,4);
    
    % Initialize min and max axes limits
    MinX = 0;MaxX = 0;
    MinY = 0;MaxY = 0;
    
    % Begin looping through saved simulation time-steps
    for TimeStepNum = 1:Sim.NumSavePoints
        % Compute coordinates of truck corner points
        p(:,1) =...
            Results.TimeStep(TimeStepNum).Truck.CG_PosVec +...
            Results.TimeStep(TimeStepNum).Truck.RotMat*...
            [Truck.CG2FrontLength Truck.Width/2]';
        p(:,2) =...
            Results.TimeStep(TimeStepNum).Truck.CG_PosVec +...
            Results.TimeStep(TimeStepNum).Truck.RotMat*...
            [Truck.CG2FrontLength -Truck.Width/2]';
        p(:,3) =...
            Results.TimeStep(TimeStepNum).Truck.CG_PosVec +...
            Results.TimeStep(TimeStepNum).Truck.RotMat*...
            [-Truck.CG2RearLength -Truck.Width/2]';
        p(:,4) =...
            Results.TimeStep(TimeStepNum).Truck.CG_PosVec +...
            Results.TimeStep(TimeStepNum).Truck.RotMat*...
            [-Truck.CG2RearLength Truck.Width/2]';
        p(:,5) = p(:,1);
        
        % Update min and max axes limits
        MinX = min(MinX,min(p(1,:)));
        MaxX = max(MaxX,max(p(1,:)));
        MinY = min(MinY,min(p(2,:)));
        MaxY = max(MaxY,max(p(2,:)));
        
        % Plot lines representing truck outer edges
        for PointNum = 1:4
            plot(...
                p(1,PointNum:PointNum + 1),...
                p(2,PointNum:PointNum + 1),...
                'k');
            hold on;
        end
        
        % Compute coordinates of trailer corner points
        p(:,1) =...
            Results.TimeStep(TimeStepNum).Truck.RearPosVec +...
            Results.TimeStep(TimeStepNum).Trailer.RotMat*...
            [0 Trailer.Width/2]';
        p(:,2) =...
            Results.TimeStep(TimeStepNum).Truck.RearPosVec +...
            Results.TimeStep(TimeStepNum).Trailer.RotMat*...
            [0 -Trailer.Width/2]';
        p(:,3) =...
            Results.TimeStep(TimeStepNum).Truck.RearPosVec +...
            Results.TimeStep(TimeStepNum).Trailer.RotMat*...
            [-Trailer.Pivot2RearLength -Trailer.Width/2]';
        p(:,4) =...
            Results.TimeStep(TimeStepNum).Truck.RearPosVec +...
            Results.TimeStep(TimeStepNum).Trailer.RotMat*...
            [-Trailer.Pivot2RearLength Trailer.Width/2]';
        p(:,5) = p(:,1);
        
        % Update min and max axes limits
        MinX = min(MinX,min(p(1,:)));
        MaxX = max(MaxX,max(p(1,:)));
        MinY = min(MinY,min(p(2,:)));
        MaxY = max(MaxY,max(p(2,:)));
        
        % Plot lines representing trailer outer edges
        for PointNum = 1:4
            plot(...
                p(1,PointNum:PointNum + 1),...
                p(2,PointNum:PointNum + 1),...
                'k');
            hold on;
        end
        hold off;
        
        % Define plot axes limits and aspect ratios
        axis([MinX MaxX MinY MaxY]);
        daspect([1 1 1]);
        pbaspect([1 1 1]);
        
        % Pause for animation purposes
        pause(Sim.SavePeriod/Sim.AnimPlayRate);
    end
end

%% Generate output plots
% Extract plot arrays from results class objects
for SampleNum = Control.NumSamples:-1:1
    PlotData.Time(SampleNum) = Results.Sample(SampleNum).Time;
    PlotData.Truck.LongVel(SampleNum) =...
        Results.Sample(SampleNum).Truck.LongVel;
    PlotData.Truck.LatVel(SampleNum) =...
        Results.Sample(SampleNum).Truck.LatVel;
    PlotData.Truck.WheelAngle(SampleNum) =...
        Results.Sample(SampleNum).Truck.WheelAngle;
    PlotData.Trailer.AngleRel2Truck(SampleNum) =...
        Results.Sample(SampleNum).Trailer.AngleRel2Truck;
end

% Plot truck forward velocity (m/s)
subplot(2,2,1);
stairs(PlotData.Time,PlotData.Truck.LongVel);
xlabel('Time (sec)');
ylabel('Truck longitudinal speed (m/s)');

% Plot truck wheel angle (deg)
subplot(2,2,2);
stairs(PlotData.Time,PlotData.Truck.WheelAngle*180/pi);
xlabel('Time (sec)');
ylabel('Truck wheel angle (deg)');

% Plot truck lateral velocity (m/s)
subplot(2,2,3);
stairs(PlotData.Time,PlotData.Truck.LatVel);
xlabel('Time (sec)');
ylabel('Truck lateral speed (m/s)');

% Plot trailer angle relative to truck (deg)
subplot(2,2,4);
stairs(PlotData.Time,PlotData.Trailer.AngleRel2Truck*180/pi);
xlabel('Time (sec)');
ylabel('Trailer angle relative to truck (deg)');