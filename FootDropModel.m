classdef FootDropModel < handle
    % Simple model of standing postural stability, consisting of foot and
    % body segments, and two muscles that create moments about the ankles,
    % tibialis anterior and soleus. 
    
    methods (Static)
        function result = soleusLength(theta)
            % Soleus length as a function of ankle angle. 
            % 
            % theta: body angle (up from prone horizontal)
            % result: soleus length
            
            origin = zeros(2,length(theta));
            for i = 1:length(theta) %this loop optionally handles a list of theta 
                origin(:,i) = [cos(theta(i)) -sin(theta(i)); sin(theta(i)) cos(theta(i))]*[.3; .03];                
            end
            insertion = [-.05; -.02]; 
            difference = origin - insertion;
            result = (difference(1,:).^2 + difference(2,:).^2).^.5;
            
            % return result with same shape as theta
            if size(theta, 1) > size(theta, 2)
                result = result';
            end
        end
        
        function result = tibialisLength(theta)
            % Soleus length as a function of ankle angle. 
            % 
            % theta: body angle (up from prone horizontal)
            % result: tibialis anterior length
            
            origin = zeros(2,length(theta));
            for i = 1:length(theta)
                origin(:,i) = [cos(theta(i)) -sin(theta(i)); sin(theta(i)) cos(theta(i))]*[.3; -.03];  
            end
            insertion = [.06; -.03];
            difference = origin - insertion;
            result = (difference(1,:).^2 + difference(2,:).^2).^.5;

            % return result with same shape as theta
            if size(theta, 1) > size(theta, 2)
                result = result';
            end
        end        
        
        function simulate(control, k, T)
            % Runs a simulation of the model and plots results. 
            % 
            % control: 0 means no control law should be usd to stabilize
            %   the model; 1 means a control law should be used
            % T: total time to simulate, in seconds
            
            restLengthS = FootDropModel.soleusLength(pi/2);
            restLengthTA = FootDropModel.tibialisLength(pi/2);
            
            S = HillTypeMuscle(16000, .6*restLengthS, .4*restLengthS);
            TA = HillTypeMuscle(2000, .6*restLengthTA, .4*restLengthTA);
            exo = Exoskeleton(k);
            
            OPTIONS = odeset('AbsTol', 1e-6, 'RelTol', 1e-5);  % use this as the final argument to ode45

            dS = .05;
            dTA = .03;
            
            theta0 = pi/2;
                        
            odefun = @(t, x) dynamics(t, x, S, TA, exo, control);            
            [time, x] = ode45(odefun, [0 T], [theta0 0 1 1], OPTIONS);
            
            fS = getForce(S, FootDropModel.soleusLength(x(:,1)), x(:,3));
%             fTA = getForce(TA, FootDropModel.tibialisLength(x(:,1)), x(:,4));
            fExo = exo.force(x(:,4));

            figure
            hold on
            subplot(2,1,1), plot(time, x(:,1))
            set(gca, 'FontSize', 18)
            ylabel('Body Angle (rad)')
            subplot(2,1,2), hold on
            plot(time, fS*dS, 'r');
            plot(time, -fExo*dTA, 'g');
            plot(time, getGravityMoment(x(:,1)), 'k')
            legend('soleus', 'tibialis+exo', 'gravity')
            set(gca, 'FontSize', 18)
            xlabel('Time (s)')
            ylabel('Torques (Nm)')
            hold off
            
        end        
    end
end

function dx_dt = dynamics(t, x, S, TA, exo, control, a_TA_in, a_S_in) 
    % Right-hand side of the dynamic equation of the model. 
    % 
    % t: time of current step (s)
    % x: model state: [ankle angle, angular velocity, soleus CE length, TA
    %   CE length]
    % S: soleus muscle object (a HillTypeModel)
    % TA: tibialis anterior muscle object (a HillTypeModel)
    
    dS = .05; % Moment arm 
    dTA = .03; % Moment arm
    
    Iankle = 90; % Moment of inertia about ankle
   
    %determine the normalized lengths of the soleus and tibialis anterior
    norm_S = S.getNormalizedLengthSE(FootDropModel.soleusLength(x(1)), x(3));
    norm_TA = TA.getNormalizedLengthSE(FootDropModel.tibialisLength(x(1)),x(4));
    
    %default activation values; used if control = 0
    a_TA = 0.4; 
    a_S = 0.05;        
        
    %control law that simulates foot drop (i.e. TA not working properly) 
    if control == 1
        if x(1) > 3*pi/4
            a_S = 0.00011;
            a_TA = 0.00000001;
        else
            a_S = 0.009;
            a_TA = 0.00000001;
        end
    end   
    
    t
    
    %calculate derivatives based on dynamic equations given
    %dx_dt is a 1D, 4 row matrix
    fexo = exo.force(norm_TA);
    x1 = x(2);
    x2 = ((S.f0M*S.forceLengthSE(norm_S)*dS) - (TA.f0M*TA.forceLengthSE(norm_TA)*dTA)-(fexo*(dTA))+ getGravityMoment(x(1)))/Iankle;
    x3 = S.getVelocity(a_S, Exoskeleton(0), x(3), norm_S);
    x4 = TA.getVelocity(a_TA, exo, x(4), norm_TA); 
    dx_dt = [x1,x2,x3,x4]';
end

function result = getGravityMoment(angle)
    % angle: angle of body segment (up from horizontal)
    % result: moment about ankle due to force of gravity on body
    
    m = 0.0145*70; % foot mass 
    lCOM = 0.5*0.25; % distance from ankle to foot centre of mass in metres
    g = 9.81; 
    result = m*g*lCOM*cos(angle-pi/2);
end

