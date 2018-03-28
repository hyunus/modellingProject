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
            
            % WRITE CODE HERE (to define odefun)
            odefun = @(t, x) HillTypeMuscle.getVelocity(afun(t), x, L-x);

            
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
        
        function result = getVelocity(a, lM, lT)
            % Calculates normalized velocity of contractile element. 
            % 
            % a: activation (between 0 and 1)
            % lM: normalized length of contractile element (this is
            %   \tilde{l}^M in the paper and the lecture slides)
            % lT: normalized length of series elastic element (this is
            %   \tilde{l}^T in the paper and the lecture slides)
            % result: normalized velocity
            
            beta = 0.1; % damping coefficient (see damped model in Millard et al.)
            
            % WRITE CODE HERE TO CALCULATE VELOCITY (use Matlab's fzero function) 
            fl = HillTypeMuscle.forceLengthCE(lM);
            fpe = HillTypeMuscle.forceLengthPE(lM);
            ft = HillTypeMuscle.forceLengthSE(lT);
            func = @(V)(a*fl*HillTypeMuscle.forceVelocityCE(V)+fpe+beta*V)-ft;
            result = fzero(func, 0);

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
            
            func = @(lT) 10*(lT - 1) + 240*(lT-1).^2;
            temp = [];            
            for i = 1:length(lT)
                if(lT(i) > 1)
                    temp(i) = func(lT(i));
                else
                    temp(i) = 0;
                end
            end            
            result = temp;           
        end
        
        function result = forceLengthPE(lM)
            % Normalized force-length curve of parallel elastic element.
            % 
            % lM: normalized length of contractile element
            % result: force produced by parallel elastic element
            
            func = @(lM) 3*(lM - 1).^2 / (0.6 + lM - 1);
            temp = [];
            for i = 1:length(lM)
                if(lM(i) > 1)
                    temp(i) = func(lM(i)); 
                else
                    temp(i) = 0;
                end
            end            
            result = temp;
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
            legend('CE', 'PE', 'SE', 'location', 'northwest')
            xlabel('Normalized length')
            ylabel('Force scale factor')
            set(gca, 'FontSize', 18)
            subplot(2,1,2)
            plot(vM, HillTypeMuscle.forceVelocityCE(vM), 'k')
            set(gca, 'FontSize', 18)
            xlabel('Normalized CE velocity')
            ylabel('Force scale factor')
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
    
xdata = [36.6602745
37.22066213
38.05621109
38.25424943
38.86712018
38.92021542
39.13190682
39.6074195
39.70315949
40.17030126
40.34367347
40.17450423
40.22484127
40.39171202
40.50548753
40.9639607
40.88305367
41.07436508
41.53385846
41.56305799
41.64829932
41.62662779
42.23920959
42.25914739
42.20453515
42.30566893
42.40680272
42.86243704
42.91536119
43.04647392
42.93269841
42.98471007
43.05405896
42.93269841
43.49230537
43.63389267
43.53950113
43.74176871
43.5915128
44.23298996
43.84290249
44.31014059
44.31620862
44.58177404
44.93110204
44.88724191
45.51919501
45.20820862
45.62285714
45.84535147
45.11718821
46.28796052
46.39147392
46.41979138
46.86680272
47.05137188
47.09363136
47.04467921
47.6581746
48.12808552
48.28684007
48.2523356
48.27256236
48.36358277
48.32312925
48.88870051
49.39177627
49.62573243
49.60752834
49.67287633
49.60752834
49.0614059
50.3087226
50.69977324
51.02651317
50.94249433
51.00991686
50.8454059
51.42793651
51.52232804
52.28959637
52.27746032
52.30780045
52.92709617
53.66991288
53.6124263
53.50209854
53.78688209
54.78131997
54.88064399
54.90019652
54.46195011
56.09020408
55.97895692
56.16099773
56.26502106
56.14582766
57.3688241
57.54653061
57.64621963
58.30503401
58.16811443
58.78758665
58.83092971
59.00285714
59.04668178
59.57426304
59.4888316
60.04453515
60.3161516
60.33580045
60.90619501
60.74714882
61.36522342
61.54557346
62.06973923
61.98630385
62.18619181
62.69171202
62.63933916
63.00438398
63.32126984
63.36794697
63.53702192
64.02087768
64.05810172
64.2921542
64.46867862
64.76025915
64.88378685
64.81804989
65.34972465
65.24955404
65.30716966
65.78893424
65.71278645
65.809161
66.11256236
66.20358277
66.63036735
66.80605118
67.14412698
67.13805896
67.15626304
67.64170522
67.67406803
67.76913379
68.41082766
68.44730807
68.66113379
68.66113379
68.73778255
69.63201814
69.71292517
69.77071591
69.80799093
69.69269841
70.32550696
70.99587949
70.96698413
70.96120505
70.96698413
71.43508908
71.29639132
71.80223023
72.4012451
72.30195011
72.70648526
72.62557823
72.69637188
72.91482086
73.33351474
73.6369161
73.88772789
74.06524743
74.01111111
74.04954195
74.06456754
75.1507292
75.03256236
75.33596372
75.34944822
75.32318892
75.57868481
76.48282086
76.36752834
76.71138322
77.7631746
78.08991453
78.87969161
];

ydata = [
    9.072908225
14.26131317
7.612785603
20.05528172
13.79733915
27.19963104
2.460015654
8.066805421
19.41754698
30.9880448
24.11307805
-1.578982849
39.99639404
11.17094875
2.550551799
49.55515812
43.40745911
19.08040219
30.58753387
24.12614052
37.78851223
9.780768344
2.739289283
53.76516467
46.81875999
42.4211036
18.30227961
28.75857342
21.66261411
38.34482526
32.63385436
13.0191494
5.97653217
-1.516735155
59.23245204
52.46774942
47.62736988
42.86272628
4.168908553
25.459686
19.8327565
37.35314286
30.67627481
64.36871934
57.72022143
48.84243188
42.43407343
26.72501636
69.42819943
54.10756816
35.23646674
62.01924023
47.22359973
75.3737894
43.70414395
80.40865623
69.19762028
53.86266727
60.73021769
46.831494
83.4796972
76.06030098
70.80445458
66.35432933
53.54338059
59.91665487
88.91260907
79.14308114
71.95994845
66.51370858
63.7717121
52.4344602
95.55337776
87.94997918
82.88077782
76.44483115
70.14800436
67.51577754
97.28153544
92.8916806
84.55703061
77.47802122
72.78412214
93.40766545
86.57457091
79.4484583
74.82510162
99.12803381
102.9344696
92.06347672
83.03028717
77.49217012
104.8121888
100.2236856
94.62512915
89.05384639
85.0017012
103.0866722
95.19209597
87.05110184
100.4145212
80.32358359
104.1993644
94.31750838
89.05854024
83.92520538
100.2303882
78.34817884
104.4312801
92.84155706
85.20714329
99.69457216
78.40176472
73.37888145
90.5055701
83.36319489
96.7004871
68.37819488
75.12801345
59.01738962
90.22659687
53.27749909
82.48370151
66.81655963
75.0841851
59.45848293
88.25359109
48.26448326
83.26077956
56.47687802
67.05923955
77.84643174
71.91143792
61.44790192
50.74535611
43.03605333
82.31694656
57.96148084
67.9147619
75.68434887
39.77553638
62.18311905
51.28379548
45.05722137
57.34639389
69.7880916
34.74931996
40.17181098
49.51396348
64.47451472
57.02436483
28.97777329
41.5102364
34.14700109
23.75927964
51.13636558
58.78448855
29.42894781
43.98684873
36.6192695
20.7843714
51.95025324
31.18988004
25.28675495
15.75425749
44.19979776
38.47656588
31.41433014
24.97566433
19.6072822
10.92056619
43.40745911
38.57136354
16.73301675
32.33626035
25.35821253
20.91777535
6.696913538
15.23644409
31.65679389
26.75419847
22.14715046
6.891246437
4.51742639
16.36744049
24.57526718
8.107835373
15.97273322
11.72031365
11.62618975
    ];

%normalize data
[M ,I] = max(ydata);
xnorm = [];
ynorm = [];

for i = 1:length(ydata)
    ynorm(i) = ydata(i)/M;
    xnorm(i) = xdata(i)/xdata(I);
end

centres = linspace(0.5, 1.5, 6);
width = .3;
lambda = .01;
result = Regression(transpose(xnorm), transpose(ynorm), centres, width, lambda, 0);
end
