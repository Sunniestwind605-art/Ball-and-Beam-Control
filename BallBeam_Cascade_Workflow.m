function results = BallBeam_Cascade_Workflow(x0, xd, doAnimation)
% BALLBEAM_CASCADE_WORKFLOW
% Cascade control for a linear ball-and-beam model.
% Outer loop: PID on ball position with a second-order setpoint prefilter
%             that cancels the outer-loop numerator (5% OS, Ts ~ 2 s target).
% Inner loop: PD on beam angle (double pole at -10).
% Order of output: simulation -> metrics -> animation -> figures.

    if nargin < 1 || isempty(x0)
        x0 = readDefaultNumber('Initial ball position x0 [m] (default 0): ', 0);
    end
    if nargin < 2 || isempty(xd)
        xd = readDefaultNumber('Desired ball position xd [m] (default -0.25): ', -0.25);
    end
    if nargin < 3 || isempty(doAnimation)
        doAnimation = true;
    end
    if ~isscalar(x0) || ~isfinite(x0) || ~isscalar(xd) || ~isfinite(xd)
        error('x0 and xd must be finite scalar positions in metres.');
    end

    %% Physical parameters and domain checks
    p.g = 9.8;
    p.L = 1.0;
    p.Mbeam = 0.24;
    p.Mball = 0.05;
    p.J0 = p.Mbeam*p.L^2/12;
    p.c = 5*p.g/7;
    p.h = 5/(7*p.g);
    p.q = p.c * p.h;   % Plant-sensor coupling factor (= 25/49)

    if abs(x0) > p.L/2 || abs(xd) > p.L/2
        error('x0 and xd must lie on the beam: between -%.3f and %.3f m.', ...
              p.L/2, p.L/2);
    end

    %% Gains (plant units: outer gains already divided by q)
    % Outer PID: target poles at -1.5 +- 1.57j and -10
    gains.outer = [34.72/p.q, 47.2/p.q, 13.0/p.q];   % [k1p, k1i, k1d]
    % Inner PD: target double pole at -10
    gains.inner = [100.0, 0.0, 20.0];                % [k2p, k2i, k2d]

    fprintf('\nBALL-AND-BEAM LINEAR CASCADE DESIGN (Full PID, Ts = 2s, OS = 5%%)\n');
    fprintf('Outer PID gains: k1p = %.2f, k1i = %.2f, k1d = %.2f\n', gains.outer);
    fprintf('Inner PD gains:  k2p = %.2f, k2i = %.2f, k2d = %.2f\n', gains.inner);

    %% Simulation of the linear plant (6 seconds runtime)
    tspan = (0:0.01:6).';
    run = simulatePlant(tspan, x0, xd, p, gains);
    run.metrics = responseMetrics(run, x0, xd);

    fprintf('\nPERFORMANCE METRICS\n');
    fprintf('Overshoot: %.2f %%\n', run.metrics.overshootPercent);
    fprintf('5%% Settling Time: %.3f s\n', run.metrics.settlingTime);
    fprintf('Peak Theta: %.3f deg\n', run.metrics.peakAngleDeg);
    fprintf('Peak Torque: %.4f Nm\n', run.metrics.peakTorqueNm);

    %% Beam animation first (figure 1)
    if doAnimation
        animateBallBeam(run, xd, p);
    end

    %% Figures appear only after the animation has finished
    plotDesignedRootLoci(p.q, gains);   % figure 10
    plotPositionResponse(run, xd);      % figure 2

    results.parameters = p;
    results.gains = gains;
    results.run = run;
end

function value = readDefaultNumber(promptText, defaultValue)
    raw = input(promptText, 's');
    if isempty(strtrim(raw))
        value = defaultValue;
    else
        value = str2double(raw);
        if ~isfinite(value)
            error('Enter a finite numeric value or press Enter for the default.');
        end
    end
end

function run = simulatePlant(tspan, x0, xd, p, gains)
    % State: [x; v; theta; omega; integralError; rf; rfDot]
    % rf is the output of the second-order setpoint prefilter.
    z0 = [x0; 0; 0; 0; 0; x0; 0];
    options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.01);
    [t,z] = ode45(@(t,z) closedLoopRhs(t,z,xd,p,gains), tspan, z0, options);
    run.time = t;
    run.x = z(:,1);
    run.velocity = z(:,2);
    run.theta = z(:,3);
    run.omega = z(:,4);
    run.integralError = z(:,5);
    run.rf = z(:,6);
    run.rfDot = z(:,7);
    [run.thetaDesired, run.motorTorque] = controllerSignals(run, xd, p, gains);
end

function dz = closedLoopRhs(~, z, xd, p, gains)
    x = z(1); v = z(2); theta = z(3); omega = z(4);
    integralError = z(5); rf = z(6); rfDot = z(7);

    % Second-order prefilter: F(s) = k1i / (k1d s^2 + k1p s + k1i)
    rfDDot = prefilterAccel(rf, rfDot, xd, gains);

    xdd = -p.c*theta;
    [thetaDesired, thetaDesiredDot] = outerLaw(x, v, xdd, integralError, rf, rfDot, rfDDot, p, gains);

    angleError = thetaDesired - theta;
    angleErrorDot = thetaDesiredDot - omega;
    angleCommand = gains.inner(1)*angleError + gains.inner(3)*angleErrorDot;
    torque = p.Mball*p.g*x + p.J0*angleCommand;

    thetadd = (torque - p.Mball*p.g*x)/p.J0;
    dz = [v; xdd; omega; thetadd; x - rf; rfDot; rfDDot];
end

function rfDDot = prefilterAccel(rf, rfDot, xd, gains)
    rfDDot = (gains.outer(2)*(xd - rf) - gains.outer(1)*rfDot) / gains.outer(3);
end

function [thetaDesired, thetaDesiredDot] = outerLaw(x, v, xdd, integralError, rf, rfDot, rfDDot, p, gains)
    e = x - rf;                 % error to the filtered setpoint
    eDot = v - rfDot;
    eDDot = xdd - rfDDot;
    thetaDesired = p.h * (gains.outer(1)*e + gains.outer(2)*integralError + gains.outer(3)*eDot);
    thetaDesiredDot = p.h * (gains.outer(1)*eDot + gains.outer(2)*e + gains.outer(3)*eDDot);
end

function [thetaDesired, torque] = controllerSignals(run, xd, p, gains)
    n = numel(run.time);
    thetaDesired = zeros(n,1);
    torque = zeros(n,1);
    for j = 1:n
        rfDDot = prefilterAccel(run.rf(j), run.rfDot(j), xd, gains);
        xdd = -p.c*run.theta(j);
        [thetaDesired(j), thetaDesiredDot] = outerLaw(run.x(j), run.velocity(j), xdd, ...
            run.integralError(j), run.rf(j), run.rfDot(j), rfDDot, p, gains);
        angleError = thetaDesired(j) - run.theta(j);
        angleErrorDot = thetaDesiredDot - run.omega(j);
        angleCommand = gains.inner(1)*angleError + gains.inner(3)*angleErrorDot;
        torque(j) = p.Mball*p.g*run.x(j) + p.J0*angleCommand;
    end
end

function metrics = responseMetrics(run, x0, xd)
    move = xd - x0;
    if abs(move) < eps
        metrics.overshootPercent = 0;
        metrics.settlingTime = 0;
    else
        progress = (run.x - x0)/move;
        metrics.overshootPercent = max(0, 100*(max(progress) - 1));
        band = 0.05*abs(move);
        lastOutside = find(abs(run.x - xd) > band, 1, 'last');
        if isempty(lastOutside)
            metrics.settlingTime = 0;
        elseif lastOutside == numel(run.time)
            metrics.settlingTime = NaN;
        else
            metrics.settlingTime = run.time(lastOutside + 1);
        end
    end
    metrics.peakAngleDeg = max(abs(run.theta))*180/pi;
    metrics.peakTorqueNm = max(abs(run.motorTorque));
end

function plotDesignedRootLoci(q, gains)
    lambda = unique([linspace(0.002,0.98,240), linspace(0.98,1.02,180), linspace(1.02,4,240)]);
    innerPoles = zeros(numel(lambda),2);
    outerPoles = zeros(numel(lambda),3);
    for k = 1:numel(lambda)
        innerPoles(k,:) = roots([1, lambda(k)*gains.inner(3), lambda(k)*gains.inner(1)]).';
        outerPoles(k,:) = roots([1, lambda(k)*q*gains.outer(3), ...
                                  lambda(k)*q*gains.outer(1), ...
                                  lambda(k)*q*gains.outer(2)]).';
    end

    figure(10); clf;
    set(gcf, 'Name', 'Cascade pole-placement root loci', 'NumberTitle', 'off', 'Color', 'w');
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot(real(innerPoles), imag(innerPoles), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(-10, 0, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 11);
    xline(0,'k:'); yline(0,'k:'); grid on; axis equal;
    xlabel('Real axis (rad/s)'); ylabel('Imaginary axis (rad/s)');
    title('Inner angle-loop root locus');

    nexttile;
    plot(real(outerPoles), imag(outerPoles), 'Color', [0 .45 .74], 'LineWidth', 1.5);
    hold on;
    plot(-1.5, 1.57, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 11);
    plot(-1.5, -1.57, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 11);
    plot(-10, 0, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 11);
    xline(0,'k:'); yline(0,'k:'); grid on; axis equal;
    xlabel('Real axis (rad/s)'); ylabel('Imaginary axis (rad/s)');
    title('Outer position PID root locus');
    sgtitle('Root Loci: Full PID Design (Ts=2s, OS=5%)');
end

function plotPositionResponse(run, xd)
    figure(2); clf;
    set(gcf, 'Name', 'Linear position response', 'NumberTitle', 'off', 'Color', 'w');
    plot(run.time, run.x, 'b-', 'LineWidth', 1.8);
    hold on;
    yline(xd, 'r--', 'Desired position', 'LineWidth', 1.3);
    grid on; xlabel('Time (s)'); ylabel('Ball position x (m)');
    title('Linear Closed-Loop Step Response (Full PID, Ts = 2s, OS = 5%)');
    legend('Actual Position x(t)', 'Desired Position', 'Location', 'best');
    hold off;
end

function animateBallBeam(run, xd, p)
    figure(1); clf;
    set(gcf, 'Color', 'w', 'Name', 'Ball and Beam Visualisation');
    for k = 1:max(1, round(length(run.time)/240)):length(run.time)
        th = run.theta(k); c = cos(th); s = sin(th);
        ends = [-p.L/2 p.L/2]; xb = ends*c; yb = ends*s;
        ballx = run.x(k)*c; bally = run.x(k)*s + 0.025;

        clf;
        plot(xb, yb, 'k-', 'LineWidth', 4); hold on;
        plot(0, 0, 'ko', 'MarkerSize', 8);
        plot(ballx, bally, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
        plot(xd*c, xd*s, 'gx', 'LineWidth', 2, 'MarkerSize', 10);
        axis equal; axis([-0.6 0.6 -0.6 0.6]); grid on;
        title(sprintf('t = %.2f s, x = %.4f m, x_d = %.4f m', run.time(k), run.x(k), xd));
        xlabel('x (m)'); ylabel('height (m)'); drawnow;
    end
end