classdef HillTypeMuscle < handle
    % Damped Hill-type muscle model adapted from Millard et al. (2013). The
    % dynamic model is defined in terms of normalized length and velocity. 
    % To model a particular muscle, scale factors are needed for force, CE
    % length, and SE length. These are given as constructor arguments. 
    
    properties (Access = public)
        f0M;               
        restingLengthCE;    
        restingLengthSE;    
    end
    
    properties (Constant)
        % curve fits for CE force-length and force-velocity 
        forceLengthRegression = getCEForceLengthRegression();
        forceVelocityRegression = getCEForceVelocityRegression();
    end
    
    methods (Access = public)
        function m = HillTypeMuscle(f0M, restingLengthCE, restingLengthSE)
            % The arguments are scale factors model is normalized 
            % f0M: maximum isometric force
            % restingLengthCE: actual length of CE (m) that corresponds to 
            %   normalized length of 1
            % restingLengthSE: % actual length of SE (m) that corresponds 
            %   to normalized length of 1
            m.f0M = f0M;
            m.restingLengthCE = restingLengthCE;
            m.restingLengthSE = restingLengthSE;
        end
        
        function result = getNormalizedLengthSE(m, muscleTendonLength, normalizedLengthCE)
            % Calculates the normalized length of the series elastic
            % element. 
            % 
            % muscleTendonLength: physical (non-normalized) length of the
            %   full muscle-tendon complex (typically found from joint 
            %   angles and musculoskeletal geometry)
            % normalizedLengthCE: normalized length of the contractile
            %   element (the state variable of the muscle model)
            % result: normalized length of the series elastic element
            
            result = (muscleTendonLength - m.restingLengthCE*normalizedLengthCE) / m.restingLengthSE;
        end
        
        function result = getForce(m, length, normalizedLengthCE)
            % length: muscle-tendon length (m)
            % normalizedLengthCE: normalized length of CE (the state
            %   variable)
            result = m.f0M * m.forceLengthSE(m.getNormalizedLengthSE(length, normalizedLengthCE));
        end
                
        function simulateIsometric(m)
            L = m.restingLengthCE + m.restingLengthSE;
            afun = @(t) t>0.5;
            %A = 1;
            odefun = @(t,x) HillTypeMuscle.getVelocity(t, x, L-x);
            
            OPTIONS = odeset('AbsTol', 1e-6, 'RelTol', 1e-5);  % use this as the final argument to ode45 
            [time, x] = ode45(odefun, [0 2], 1, OPTIONS);
            
            figure
            subplot(2,1,1)
            plot(time, x*m.restingLengthCE)
            ylabel('CE length')
            set(gca, 'FontSize', 18)
            subplot(2,1,2)
            plot(time, m.f0M*HillTypeMuscle.forceLengthSE(m.getNormalizedLengthSE(L, x)));
            xlabel('Time (s)')
            ylabel('Force')
            set(gca, 'FontSize', 18)
        end       
    end
    
    methods (Static)
        
        function result = getVelocity(a,exo,lM,lT)
            % Calculates normalized velocity of contractile element. 
            % 
            % a: activation (between 0 and 1)
            % lM: normalized length of contractile element (this is
            %   \tilde{l}^M in the paper and the lecture slides)
            % lT: normalized length of series elastic element (this is
            %   \tilde{l}^T in the paper and the lecture slides)
            % result: normalized velocity
            
%             if (t=<.5)
              %  a=0;
            %end
            %if (t>0.5)
%                 a=t;
%             end
     
            beta = 0.1; % damping coefficient (see damped model in Millard et al.)
            % WRITE CODE HERE TO CALCULATE VELOCITY (use Matlab's fzero function) 
            fT = HillTypeMuscle.forceLengthSE(lT);
            fPE = HillTypeMuscle.forceLengthPE(lM);
            fL = HillTypeMuscle.forceLengthCE(lM);
            fEXO = exo.force(lM);
            func = @(vM) (1*(((a*fL*HillTypeMuscle.forceVelocityCE(vM))+ (fPE) + (beta*vM)+ fEXO)*cos(0)) - (1*fT));
            vMNorm = fzero(func, 0);
            result = vMNorm;
        end
        
        %for plotting purposes
        function result = forceExo(k,lM)
             for i = 1:length(lM)
                if(lM(i) < 1)
                    result(i) = 0;    
                else
                    result(i) = k*(lM(i)-1);
                end
            end
        
        end
        
        function result = forceLengthCE(lM)
            % Normalized force-length curve of contractile element. 
            % 
            % lM: contracile element length
            % result: force-length scale factor
            
            result = HillTypeMuscle.forceLengthRegression.eval(lM);            
        end
        
        function result = forceVelocityCE(vM)
            % Normalized force-velocity curve of contractile element.  
            % 
            % vM: contracile element velocity
            % result: force-velocity scale factor
            
            result = max(0, HillTypeMuscle.forceVelocityRegression.eval(vM));
        end
        
        function result = forceLengthSE(lT)
            % Normalized force-length curve of series elastic element. 
            % 
            % lT: normalized length of tendon (series elastic element)
            % result: force produced by tendon
            
            % WRITE CODE HERE
            for i = 1:length(lT)
                if(lT(i) < 1)
                    result(i) = 0;    
                else
                    result(i) = 10*(lT(i)-1)+240*(lT(i)-1).^2;
                end
            end
        end
        
        function result = forceLengthPE(lM)
            % Normalized force-length curve of parallel elastic element.
            % 
            % lM: normalized length of contractile element
            % result: force produced by parallel elastic element
            
            % WRITE CODE HERE   
            
              for i = 1:length(lM)
                if(lM(i) < 1)
                    result(i) = 0;    
                else
                    result(i) = (3*(lM(i)-1).^2)/(0.6+lM(i)-1);
                end
            end
        end

        
        function plotCurves()
            % Plot force-length, force-velocity, SE, and PE curves. 
            
            lM = 0:.01:1.8;
            vM = -1.2:.01:1.2;
            lT = 0:.01:1.07;
            figure
            subplot(2,1,1), hold on
            plot(lM, HillTypeMuscle.forceLengthCE(lM), 'r')
            plot(lM, HillTypeMuscle.forceLengthPE(lM), 'g')
            plot(lT, HillTypeMuscle.forceLengthSE(lT), 'b')
            plot(lM, HillTypeMuscle.forceExo(0.4,lM),'y')
            legend('CE', 'PE', 'SE', 'location', 'northwest')
            xlabel('Normalized length')
            ylabel('Force scale factor')
            set(gca, 'FontSize', 18)
            axis([0 2 0 1.1])  
            subplot(2,1,2)
            plot(vM, HillTypeMuscle.forceVelocityCE(vM), 'k')
            set(gca, 'FontSize', 18)  
            xlabel('Normalized CE velocity')
            ylabel('Force scale factor')
            axis([-1.5 1.5 0 1.5]) 
        end
        
    end
    
end

function result = getCEForceVelocityRegression()
    % result: regression model of contractile element force-length curve
    %   based on data in Millard et al.
    
    data = [-1.0028395556708567, 0.0024834319945283845
    -0.8858611825192801, 0.03218792009622429
    -0.5176245843258415, 0.15771090304473967
    -0.5232565269687035, 0.16930496922242444
    -0.29749770052593094, 0.2899790099290114
    -0.2828848376217543, 0.3545364496120378
    -0.1801231103040022, 0.3892195938775034
    -0.08494610976156225, 0.5927831890757294
    -0.10185137142991896, 0.6259097662790973
    -0.0326643239546236, 0.7682365981934388
    -0.020787245583830716, 0.8526638522676352
    0.0028442725407418212, 0.9999952831301149
    0.014617579774061973, 1.0662107025777694
    0.04058866536166583, 1.124136223202283
    0.026390887007381902, 1.132426122025424
    0.021070257776939272, 1.1986556920827338
    0.05844673474682183, 1.2582274002971627
    0.09900238201929201, 1.3757434966156459
    0.1020023112662436, 1.4022310794556732
    0.10055894908138963, 1.1489210160137733
    0.1946227683309354, 1.1571212943090965
    0.3313459588217258, 1.152041225442796
    0.5510200231126625, 1.204839508502158];

    velocity = data(:,1);
    force = data(:,2);

    centres = linspace(-.5, 0, 6);
    width = .3;
    lambda = .01;
    result = Regression(velocity, force, centres, width, lambda, 1);
end

function result = getCEForceLengthRegression()
    % result: Regression model of normalized CE force-length curve. Data 
    %   from Winters et al. (2011) Fig. 3C, normalized so that max force is
    %   ~1 and length at max force is ~1. 
    
    % WRITE CODE HERE (use WebPlotDigitizer to extract force-length points 
    % from Winters et al. (2011) Figure 3C. Click "View Data", select all, 
    % cut, and paste below, then write code to normalize and perform the 
    % regression curve fit. Use the plotCurves() function to verify that 
    % the result is reasonable. 
    % 
    % Don't forget to normalize data so optimal length = 1 and peak = 1.
    
    data = [ 
            37.34089116719243, 9.90000000000002
            38.36928233438486, 14.700000000000017
            39.24329652996846, 24.60000000000001
            39.38742113564669, 3.9000000000000057
            40.34286277602524, 17.700000000000017
            40.42448738170347, 21.90000000000002
            40.43414826498423, 36.60000000000001
            41.36612776025237, 14.700000000000017
            41.36671924290221, 15.600000000000023
            41.37736593059937, 31.800000000000026
            41.453075709779185, 27.000000000000014
            41.75216876971609, 2.1000000000000227
            42.0871451104101, 31.800000000000026
            42.48836750788644, 42.30000000000001
            42.80678233438486, 46.800000000000026
            42.807768138801265, 48.30000000000001
            42.87066246056783, 24.00000000000003
            43.124605678233436, 50.40000000000002
            43.35745268138801, 44.70000000000002
            43.42251577287067, 23.700000000000017
            43.430007886435334, 35.10000000000002
            43.430007886435334, 35.10000000000002
            43.442429022082024, 54.000000000000014
            43.76005520504732, 57.30000000000001
            43.89471608832808, 22.20000000000003
            44.383477917981075, 45.90000000000002
            44.47200315457414, 60.60000000000001
            45.25650630914827, 54.30000000000001
            45.41384069400631, 53.70000000000002
            45.565063091482656, 43.800000000000026
            45.56683753943218, 46.500000000000014
            45.73856466876972, 67.80000000000001
            45.97712933753944, 70.80000000000001
            46.444992113564666, 62.70000000000002
            46.45208990536278, 73.50000000000001
            46.52977129337539, 71.70000000000002
            46.748619873817034, 44.70000000000002
            46.768927444794954, 75.60000000000002
            47.16660094637224, 80.70000000000002
            47.47022870662461, 62.70000000000002
            47.476143533123036, 71.70000000000002
            47.640575709779185, 81.9
            47.70958201892745, 66.9
            48.193611987381715, 83.4
            48.98067823343849, 81
            48.98146687697161, 82.20000000000002
            48.9836356466877, 85.50000000000001
            49.04771293375395, 63.000000000000014
            49.529574132492115, 76.20000000000002
            49.69440063091483, 87
            49.769913249211356, 81.9
            49.85094637223975, 85.20000000000002
            50.32570977917982, 87.60000000000001
            50.63249211356468, 74.4
            50.71845425867508, 85.20000000000002
            50.719834384858046, 87.30000000000001
            50.72200315457414, 90.60000000000001
            50.792586750788644, 78
            51.03055993690852, 80.10000000000001
            51.1875, 78.9
            51.51005520504732, 89.70000000000002
            51.74743690851736, 90.9
            52.69262618296531, 89.10000000000001
            53.32334384858044, 88.80000000000001
            53.553430599369094, 78.9
            53.5563880126183, 83.4
            53.56368296529969, 94.5
            53.64392744479496, 96.60000000000001
            53.72003154574133, 92.4
            53.88071766561515, 96.9
            53.956427444794954, 92.10000000000001
            54.194006309148264, 93.60000000000001
            54.35232649842272, 94.5
            54.75000000000001, 99.60000000000001
            55.77287066246057, 96
            56.16719242902209, 96
            56.48501577287066, 99.60000000000001
            56.87933753943218, 99.60000000000001
            57.19479495268139, 99.60000000000001
            57.272673501577295, 98.10000000000001
            57.5891167192429, 99.60000000000001
            57.89905362776026, 91.20000000000002
            57.90437697160884, 99.30000000000001
            58.45070977917982, 90.60000000000001
            58.61218454258676, 96.30000000000001
            58.929416403785496, 99
            59.39747634069401, 91.20000000000002
            59.401813880126184, 97.80000000000001
            59.71589116719243, 95.70000000000002
            59.87440851735016, 96.9
            60.031348580441644, 95.70000000000002
            60.660883280757105, 93.60000000000001
            60.66482649842272, 99.60000000000001
            61.28608044164039, 84.9
            61.369479495268145, 91.80000000000001
            61.371451104100956, 94.80000000000001
            61.43848580441641, 76.80000000000001
            61.44026025236593, 79.50000000000001
            61.44578075709779, 87.9
            62.23541009463723, 89.4
            62.397673501577295, 96.30000000000001
            62.62342271293376, 79.80000000000001
            63.17981072555206, 86.4
            63.3945189274448, 53.10000000000001
            63.412657728706634, 80.70000000000002
            63.477720820189276, 59.70000000000002
            63.490733438485805, 79.50000000000001
            63.65240536277603, 85.50000000000001
            63.734029968454266, 89.70000000000002
            63.81111987381704, 87
            63.88564668769717, 80.4
            63.96175078864354, 76.20000000000002
            64.5964116719243, 81.9
            64.59976340694007, 87
            64.81348580441642, 52.20000000000002
            64.8144716088328, 53.70000000000002
            65.14471608832808, 76.20000000000002
            65.37342271293376, 64.20000000000002
            65.61573343848582, 72.9
            65.75709779179812, 48.000000000000014
            65.76932176656152, 66.60000000000001
            65.77030757097793, 68.10000000000001
            66.16738958990537, 72.30000000000001
            66.48481861198738, 75.30000000000001
            66.95208990536278, 66.30000000000001
            67.01577287066246, 43.20000000000002
            67.26557570977918, 63.30000000000001
            67.40536277602524, 36.00000000000003
            67.57334384858045, 51.60000000000002
            67.89609621451106, 62.70000000000002
            68.44617507886436, 59.70000000000002
            68.58300473186121, 27.900000000000006
            69.45958201892745, 41.70000000000002
            70.24053627760253, 30.00000000000003
            70.31940063091483, 30.00000000000003
            70.64708201892745, 48.60000000000002
            71.42665615141956, 34.80000000000001
            72.44538643533124, 24.90000000000002
            72.9191640378549, 25.800000000000026
            73.39826498422714, 34.80000000000001
            73.46293375394322, 13.200000000000017
            73.46628548895899, 18.30000000000001
            73.46687697160883, 19.200000000000017
            74.48817034700316, 13.200000000000017
            75.28016561514197, 18.30000000000001
            75.51320977917982, 12.90000000000002
            76.45682176656152, 8.700000000000017
    ];
    length = data(:,1);
    force = data(:,2);
    [maxforce,i] = max(force);
    normlength = length/length(i);
    normforce = force/maxforce;
    centres = linspace(normlength(1), normlength(end), 50);
    width = 0.158;
    lambda = .01;
    result = Regression(normlength, normforce, centres, width, lambda, 1);
end
