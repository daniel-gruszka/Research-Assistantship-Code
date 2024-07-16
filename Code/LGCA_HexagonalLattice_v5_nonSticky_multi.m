function LGCA_HexagonalLattice_v5_nonSticky_multi(JFunction,J0Function,KFunction,nameOfVideo,VariableDataName)
    %{ JFunction- Value of J
    %J0Function - Value of J0
    %KFunction - Value of K
    %nameOfVideo - Name of the video being produced
    %VariableDataName - Name for file that variable data will be stored in
    %OUTPUT FILES:
    %1 MP4 video of simulation that is named based on variable nameOfVideo
    %296 matlab data files called nameOfVide0_(J/J0)_DataImage(iteration #).
    %   The J/J0 is dependent on whether it is the J or J0 simulation. So
    %   there is 148 of each corresponding to the first 50 iteration and then
    %   every 20 after that. These are the images that generate aggregate
    %   count, radius of gyration, circularity and standard deviation of
    %   the area.
    %296 matlab data files called
    %   nameOfVideo_(J/J0)_uncleanedAggCount(iteration #). This is the same
    %   as previous saved files but is only used for the uncleaned
    %   aggregate counting. 
    %592 matlab data files called
    %   nameOfVideo(Particles/Clocks)(J/J0)(iteration #). These are the saved
    %   particles or clock data for J or J0. These are saved for the first
    %   50 iterations and then every 20 after. 
    %2 matlab data files called VariableDataName(J/J0). This is the file
    %   that has saved all the data for either the J or J0 simulation. The
    %   columns are as follows:
    %   r | r Local | Mean Square Data | Aggregate Count based on cleaned
    %   data | Aggregate Count based on Unclean data | Radius of Gyration |
    %   Average Circularity | Standard Deviation of Area
    clc
    %clf
    tic
    close all
    %Grid size
    xNodes = 50;
    yNodes = 50;
    %Percentage of possible particle locations filled
    percentCapacity = 0.4;
    %Plots every 10 iterations 
    % % % plotTen = false;
    %Calculates the Mean Square Displacement of the Particles
    msdTest = true;
    %Calculates R and RLocal
    radiusTest = true;
    %Particles have a sticky phase
    sticky = false; 
    %Wait Time between distribution
    waitTime = 0;
    %Dot Size for main simulation
    dotSize = 6;
    %The speed of the clocks (how much they baseline change each iteration)
    clockSpeed = 0;
    %How delayed the clocks are from synchronizing with neighbors
    clockDelay =0;
    %Number of Iterations
    numIterations = 2000;

    
    %How many times clocks update per LGCA step
    clockStep = 10;
    
    J=JFunction;
    J0=J0Function;
    K=KFunction;

    %Names for saving data and movie
    saveParticles = [nameOfVideo, 'Particles'];
    saveClocks = [nameOfVideo, 'Clocks'];
    videoName = nameOfVideo;

    
    

    %Initialize particle array
    particlesJ = zeros(xNodes,yNodes,6);
    
    
    %Keeping track of where the particles start for MSD
    xInitialJ = zeros(xNodes,yNodes,6);
    yInitialJ = zeros(xNodes,yNodes,6);
    
    %Adding the clock
    thetaJ = zeros(xNodes,yNodes,6);
   
    %Initializing msd, r, rlocal, aggCount, uncleaned agg count, radius of
    %gyration and std of area
    msdDataYJ = zeros(numIterations,1);
    msdDataYJ0=msdDataYJ;
    rDataYJ = zeros(numIterations,1);
    rDataYJ0=rDataYJ;
    rLocalDataJ= zeros(numIterations,1);
    rLocalDataJ0=rLocalDataJ;
    aggCountYJ = zeros(148,1);
    aggCountYJ0 = aggCountYJ;
    aggCountUnCleanYJ = zeros(148,1);
    aggCountUnCleanYJ0 = zeros(148,1);
    radiusOfGJ= zeros(148,1);
    radiusOfGJ0 = radiusOfGJ;
    avgCircJ = zeros(148,1);
    avgCircJ0 = avgCircJ;
    stdAreaJ = zeros(148,1);
    stdAreaJ0 = stdAreaJ;
    

    v = VideoWriter(videoName,'MPEG-4');
    open(v);
    
    
    numParticles = 0;    
    %sets random points as particles
    while numParticles < xNodes*yNodes*6*percentCapacity
        iInitial = randi([1,xNodes]);
        jInitial = randi([1,yNodes]);
        kInitial = randi([1,6]);
        if(particlesJ(iInitial,jInitial,kInitial)==0)
            particlesJ(iInitial,jInitial,kInitial) =1;
            thetaJ(iInitial,jInitial,kInitial)=2*pi*rand+50*pi;
            numParticles = numParticles +1;
        end
    end
    
    
    %Calculate initial particle location for MSD
    for iMSD =1:xNodes
        for jMSD = 1:yNodes
            for kMSD = 1:6
               
                if(particlesJ(iMSD,jMSD,kMSD)==1)
                    xInitialJ(iMSD,jMSD,kMSD) = iMSD + jMSD*cos(pi/3);
                    yInitialJ(iMSD,jMSD,kMSD) = jMSD*sin(pi/3);
                    
                end
            end
        end
    end
    %Copy everything for J0 so J and J0 have same starting conditions.
    particlesJ0=particlesJ;
    xInitialJ0=xInitialJ;
    yInitialJ0=yInitialJ;
    thetaJ0=thetaJ;

    
    %Plot initial particles
    f1 = figure;
    f2 = figure;
    plotParticles2(particlesJ,xNodes,yNodes,thetaJ,dotSize,numParticles,true);
    plotParticles2(particlesJ0,xNodes,yNodes,thetaJ0,dotSize,numParticles,false);
    
   
    
    
    figure(f2)
    cleanPar = cleanAggregates(xNodes,yNodes,particlesJ0); %Clean Aggregates
    plotParticlesSimple(cleanPar,xNodes,yNodes,30); %Plots Black and White cleaned aggregates
    F = getframe(gcf);
    X = frame2im(F);
    X1 = X(:,:,1);
    save([nameOfVideo,'_J_DataImage',num2str(1),'.mat'],"X1") %save B/W Image
    aggCountYJ(1)=countAggregates(X1); %Counts based on cleaned aggregates
    [radiusOfGJ(1), avgCircJ(1),stdAreaJ(1)] = radiusOfG(X1); %Calculates radius of Gyration, Circularity and std area based on cleaned aggregates

    %Same as previous but now for J0 instead of J
    cleanPar = cleanAggregates(xNodes,yNodes,particlesJ0);
    plotParticlesSimple(cleanPar,xNodes,yNodes,30);
    F = getframe(gcf);
    X = frame2im(F);
    X1 = X(:,:,1);
    save([nameOfVideo,'_J0_DataImage',num2str(1),'.mat'],"X1")
    aggCountYJ0(1)=countAggregates(X1);
    [radiusOfGJ0(1), avgCircJ0(1),stdAreaJ0(1)] = radiusOfG(X1); 

    %Now counts based on uncleaned data and saves that image
    plotParticlesSimple(particlesJ,xNodes,yNodes,30);
    F = getframe(gcf);
    X = frame2im(F);
    X1 = X(:,:,1);
    save([nameOfVideo,'_J_uncleanedAggCount',num2str(1),'.mat'],"X1")
    aggCountUnCleanYJ(1)=countAggregates(X1);
    
    %Same but for J0 instead of J
    plotParticlesSimple(particlesJ0,xNodes,yNodes,30);
    F = getframe(gcf);
    X = frame2im(F);
    X1 = X(:,:,1);
    save([nameOfVideo,'_J0_uncleanedAggCount',num2str(1),'.mat'],"X1")
    aggCountUnCleanYJ0(1)=countAggregates(X1);

    %Save initial data
    save([saveParticles, 'J', num2str(1),'.mat'], "particlesJ")
    save([saveParticles, 'J0', num2str(1),'.mat'], "particlesJ0")
    save([saveClocks, 'J', num2str(1),'.mat'], "thetaJ")
    save([saveClocks, 'J0', num2str(1),'.mat'], "thetaJ0")
    
    %calculate first r and rlocal for J and J0
    rMagJ = radius(particlesJ,xNodes,yNodes,thetaJ);
    localRJ = localRvalue(thetaJ,xNodes,yNodes);
    rDataYJ(1) = rMagJ;
    rLocalDataJ(1)=localRJ;

    rMagJ0 = radius(particlesJ0,xNodes,yNodes,thetaJ0);
    localRJ0 = localRvalue(thetaJ0,xNodes,yNodes);
    rDataYJ0(1) = rMagJ0;
    rLocalDataJ0(1)=localRJ0;

    %Calculate MSD for J and J0
    msdJ = MSD(particlesJ,xNodes,yNodes,xInitialJ,yInitialJ);
    msdDataYJ(1) = msdJ;
    msdJ0 = MSD(particlesJ0,xNodes,yNodes,xInitialJ0,yInitialJ0);
    msdDataYJ0(1) = msdJ0;

    figure(f1)
    %Annotation for the figure
    annotation('textbox', [0.1, 0.87, 0.1, 0.1], 'String', "Iteration: " + 1)
    annotation('textbox', [0.2, 0.87, 0.1, 0.1], 'String', "Time: " +toc)

    annotation('textbox', [0.1, 0.81, 0.1, 0.1], 'String', "J= " +J)
    annotation('textbox', [0.3, 0.87, 0.1, 0.1], 'String', "K= " +K)
    annotation('textbox', [0.1,0.75, 0.1, 0.1], 'String', "r= " + rMagJ)
    annotation('textbox', [0.2, 0.75, 0.1, 0.1], 'String', "localR= " + localRJ)
    annotation('textbox', [0.3, 0.75, 0.1, 0.1], 'String', "Count= " + aggCountYJ(1))
    annotation('textbox', [0.4, 0.75, 0.1, 0.1], 'String', "nr^2= " + aggCountYJ(1)*rMagJ^2)
    
    annotation('textbox', [0.55, 0.81, 0.1, 0.1], 'String', "J0= " + J0)
    annotation('textbox', [0.55,0.75, 0.1, 0.1], 'String', "r= " + rMagJ0)
    annotation('textbox', [0.65, 0.75, 0.1, 0.1], 'String', "localR= " + localRJ0)
    annotation('textbox', [0.75, 0.75, 0.1, 0.1], 'String', "Count= " + aggCountYJ0(1))
    annotation('textbox', [0.85, 0.75, 0.1, 0.1], 'String', "nr^2= " + aggCountYJ0(1)*rMagJ0^2)
    
    %Write to the video
    frame = getframe(1); 
    writeVideo(v,frame);

    dataCounter = 2;
    for ii = 2:numIterations
        
        %Moves Particles
        [particlesJ,xInitialJ,yInitialJ,thetaJ] = moveParticles(particlesJ,xNodes,yNodes,xInitialJ,yInitialJ,thetaJ);
        [particlesJ0,xInitialJ0,yInitialJ0,thetaJ0] = moveParticles(particlesJ0,xNodes,yNodes,xInitialJ0,yInitialJ0,thetaJ0);
        
        %Calculate MSD and save data for later 
        if(msdTest)
            msdJ = MSD(particlesJ,xNodes,yNodes,xInitialJ,yInitialJ);
            msdDataYJ(ii) = msdJ;
            msdJ0 = MSD(particlesJ0,xNodes,yNodes,xInitialJ0,yInitialJ0);
            msdDataYJ0(ii) = msdJ0;
        end
        %Updates clocks
        for iClockStep = 1:clockStep
            thetaJ = kuramoto(particlesJ,xNodes,yNodes,thetaJ,K/clockStep,clockSpeed/clockStep,clockDelay);
            thetaJ0 = kuramoto(particlesJ0,xNodes,yNodes,thetaJ0,K/clockStep,clockSpeed/clockStep,clockDelay);
        end
        
        %interaction step
        [particlesJ,thetaJ]= ALGCAClocks2(particlesJ,xNodes,yNodes,thetaJ,J,0,sticky);
        [particlesJ0,thetaJ0]= ALGCAClocks2(particlesJ0,xNodes,yNodes,thetaJ0,0,J0,sticky);
        
        %%%[particles,theta,kTracker]= SLGCAClocks(particles,xNodes,yNodes,theta,J,J0,xTracker,yTracker,kTracker,sticky);
        
        if(radiusTest)
            rMagJ = radius(particlesJ,xNodes,yNodes,thetaJ);
            localRJ = localRvalue(thetaJ,xNodes,yNodes);
            rDataYJ(ii) = rMagJ;
            rLocalDataJ(ii)=localRJ;

            rMagJ0 = radius(particlesJ0,xNodes,yNodes,thetaJ0);
            localRJ0 = localRvalue(thetaJ0,xNodes,yNodes);
            rDataYJ0(ii) = rMagJ0;
            rLocalDataJ0(ii)=localRJ0;
        end
        
        
        %Plots Moved Particles
        clf
        hold on
        plotParticles2(particlesJ,xNodes,yNodes,thetaJ,dotSize,numParticles,true);
        plotParticles2(particlesJ0,xNodes,yNodes,thetaJ0,dotSize,numParticles,false);

        %Calculates count, radius of gyration, circularity, std area for
        %cleaned aggregates and count for uncleaned aggregates for the
        %first 50 iterations and then every 20 after that
        if(ii<=50 || (0==mod(ii,20)))

            %Count, radius of gyration, circularity and std area for 
            % cleaned aggregates J    
            cleanParJ = cleanAggregates(xNodes,yNodes,particlesJ);
            plotParticlesSimple(cleanParJ,xNodes,yNodes,30);
            F = getframe(gcf);
            X = frame2im(F);
            X1 = X(:,:,1);
            save([nameOfVideo,'_J_DataImage',num2str(ii),'.mat'],"X1")
            aggCountYJ(dataCounter)=countAggregates(X1);
            [radiusOfGJ(dataCounter), avgCircJ(dataCounter),stdAreaJ(dataCounter)] = radiusOfG(X1);
    
            %Count, radius of gyration, circularity and std area for
            %cleaned aggregates J0
            cleanPar = cleanAggregates(xNodes,yNodes,particlesJ0);
            plotParticlesSimple(cleanPar,xNodes,yNodes,30);
            F = getframe(gcf);
            X = frame2im(F);
            X1 = X(:,:,1);
            save([nameOfVideo,'_J0_DataImage',num2str(ii),'.mat'],"X1")
            aggCountYJ0(dataCounter)=countAggregates(X1);
            [radiusOfGJ0(dataCounter), avgCircJ0(dataCounter),stdAreaJ0(dataCounter)] = radiusOfG(X1);

            %Calculates count for uncleaned aggregates J
            plotParticlesSimple(particlesJ,xNodes,yNodes,30);
            F = getframe(gcf);
            X = frame2im(F);
            X1 = X(:,:,1);
            save([nameOfVideo,'_J_uncleanedAggCount',num2str(ii),'.mat'],"X1")
            aggCountUnCleanYJ(dataCounter)=countAggregates(X1);

            %Calculates count for uncleaned aggregates J0
            plotParticlesSimple(particlesJ0,xNodes,yNodes,30);
            F = getframe(gcf);
            X = frame2im(F);
            X1 = X(:,:,1);
            save([nameOfVideo,'_J0_uncleanedAggCount',num2str(ii),'.mat'],"X1")
            aggCountUnCleanYJ0(dataCounter)=countAggregates(X1);


            %Saves particle and clock data
            save([saveParticles, 'J', num2str(ii),'.mat'], "particlesJ")
            save([saveParticles, 'J0', num2str(ii),'.mat'], "particlesJ0")
            save([saveClocks, 'J', num2str(ii),'.mat'], "thetaJ")
            save([saveClocks, 'J0', num2str(ii),'.mat'], "thetaJ0")
            dataCounter = dataCounter+1; 
        
        end
        %Annotation for the figure
        figure(f1)
        annotation('textbox', [0.1, 0.87, 0.1, 0.1], 'String', "Iteration: " + ii)
        annotation('textbox', [0.2, 0.87, 0.1, 0.1], 'String', "Time: " +toc)

        annotation('textbox', [0.1, 0.81, 0.1, 0.1], 'String', "J= " +J)
        annotation('textbox', [0.3, 0.87, 0.1, 0.1], 'String', "K= " +K)
        annotation('textbox', [0.1,0.75, 0.1, 0.1], 'String', "r= " + rMagJ)
        annotation('textbox', [0.2, 0.75, 0.1, 0.1], 'String', "localR= " + localRJ)
        annotation('textbox', [0.3, 0.75, 0.1, 0.1], 'String', "Count= " + aggCountYJ(dataCounter-1))
        annotation('textbox', [0.4, 0.75, 0.1, 0.1], 'String', "nr^2= " + aggCountYJ(dataCounter-1)*rMagJ^2)

        annotation('textbox', [0.55, 0.81, 0.1, 0.1], 'String', "J0= " + J0)
        annotation('textbox', [0.55,0.75, 0.1, 0.1], 'String', "r= " + rMagJ0)
        annotation('textbox', [0.65, 0.75, 0.1, 0.1], 'String', "localR= " + localRJ0)
        annotation('textbox', [0.75, 0.75, 0.1, 0.1], 'String', "Count= " + aggCountYJ0(dataCounter-1))
        annotation('textbox', [0.85, 0.75, 0.1, 0.1], 'String', "nr^2= " + aggCountYJ0(dataCounter-1)*rMagJ0^2)
        
        %Write to video
        frame = getframe(1); 
        writeVideo(v,frame); 
        pause(waitTime)
    end
    %Closes video
    close(v) 
    
    close
    %Prepares and saves data for J
    %Final data has columns: r, rLocal, MSD, Agg Count Cleaned, Agg Count
    %Uncleaned, Radius of Gyration, Circularity, STD Area
    aggCountJ =[aggCountYJ;zeros(numIterations-length(aggCountYJ),1) ];
    aggCountUnCleanYJ = [aggCountUnCleanYJ; zeros(numIterations-length(aggCountUnCleanYJ),1) ];
    radiusOfGJ = [radiusOfGJ;zeros(numIterations-length(aggCountYJ),1)  ];
    avgCircJ = [avgCircJ;zeros(numIterations-length(aggCountYJ),1)  ];
    stdAreaJ = [stdAreaJ;zeros(numIterations-length(stdAreaJ),1)  ];
    savedDataJ = [rDataYJ rLocalDataJ msdDataYJ aggCountJ aggCountUnCleanYJ radiusOfGJ avgCircJ stdAreaJ];
    save([VariableDataName,'_J.mat'],"savedDataJ")

    %Perpares and saves data for J0
    aggCountJ0 =[aggCountYJ0;zeros(numIterations-length(aggCountYJ0),1) ];
    aggCountUnCleanYJ0 = [aggCountUnCleanYJ0; zeros(numIterations-length(aggCountUnCleanYJ0),1) ];
    radiusOfGJ0 = [radiusOfGJ0;zeros(numIterations-length(aggCountYJ),1)  ];
    avgCircJ0 = [avgCircJ0;zeros(numIterations-length(aggCountYJ),1)  ];
    stdAreaJ0 = [stdAreaJ0;zeros(numIterations-length(stdAreaJ0),1)  ];
    savedDataJ0 = [rDataYJ0 rLocalDataJ0 msdDataYJ0 aggCountJ0 aggCountUnCleanYJ0 radiusOfGJ0 avgCircJ0 stdAreaJ0];
    save([VariableDataName,'_J0.mat'],"savedDataJ0")

    %Saves final figure
    saveFigName = [num2str(J),'J_',num2str(J0),'J0_',num2str(K),'K_multi.fig'];
    savefig(saveFigName);


    
    
    
    %Functions
    
    function plotParticles2(p,x,y,theta,dotSize,numParticles,JorJ0)
        %Main Plot function for displaying the data
        %p - particle data (3 dimension array (xNodes by yNodes by 6))
        %x - x grid size (Int)
        %y - y grid size (Int)
        %theta- clock data (3 dimension array (xNodes by yNodes by 6))
        %dotSize - Size of dots to be ploted (Int)
        %numParticles - number of particles in simulation (Int)
        %JorJ0 - tells whether we are ploting J or J0 data (Boolean)

        stencil = [(1/3),0; 0,(1/3);-(1/3),(1/3);-(1/3),0;0,-(1/3);(1/3),-(1/3);0,0;]; %Used for off setting points in hexagons
        %initializing the points and RGB
        Points = zeros(numParticles,2);
        RGB = zeros(numParticles,3);
        pointCounter = 1;
        for i =1:x
            for j = 1:y
                for k = 1:6
                    if p(i,j,k)==1
                        
                        %Calculating where each particle needs to be placed and
                        %its RGB color
                        temp = stencil(k,:);
                        if(i>y-floor(j/2)) %Changing the shape of the layout to a rectangle for display purposes 

                            %Calculating the cartesian coordinates of the
                            %points
                            Points(pointCounter,1)= (i-x)+temp(1)+ (j+temp(2))*cos(pi/3);
                            Points(pointCounter,2)= (j+temp(2))*sin(pi/3);
                            
                        else

                            %Calculating the cartesian coordinates of the
                            %points
                            Points(pointCounter,1)= i+temp(1)+ (j+temp(2))*cos(pi/3);
                            Points(pointCounter,2)= (j+temp(2))*sin(pi/3);
                            
                        end
                        
                        %Calculates the color of point based on clock
                        RGB(pointCounter,:) = hsv2rgb(mod(theta(i,j,k),2*pi)/(2*pi),1,1);
                        pointCounter = pointCounter + 1;
        
                    end
                end
            end
        end
        
        figure(f1)

        axis off
        if(JorJ0)
            clf
            subplot(1,2,1)
        else
            subplot(1,2,2)
        end

        %PLot points
        scatter(Points(:,1),Points(:,2),dotSize,RGB,'filled');
        
        f = gcf;
        f.Position = [0 0 1400 775];
        set(gcf,'color','w');
        xlim([0 x+1])
        pbaspect([1.1,1,1])
        ylim([0 y])
        set(gca,'LooseInset',get(gca,'TightInset'));
        axis off
    end
    
    function msd = MSD(p,x,y,xInt,yInt)
        %Calculates the Mean Square Displacement of the particles
        %p - particle data (3 dimension array (xNodes by yNodes by 6))
        %x - x size of gird (Int)
        %y - y size of grid (Int) 
        %xInt - the inital x coordinates of the particles (3 dimension array (xNodes by yNodes by 6))
        %yInt - the inital y coordinates of the particles (3 dimension array (xNodes by yNodes by 6))
        difSum=0;
        counter =0;
        for i = 1:x
            for j = 1:y
                for k = 1:6
                    if(p(i,j,k)==1)
                        %Calculates the displacement and sums it
                        displacement = sqrt((i + j*cos(pi/3)-xInt(i,j,k))^2 + (j*sin(pi/3) - yInt(i,j,k))^2);
                        difSum = difSum + displacement;
                        counter = counter +1;
                    end
                end
            end
        end
        msd = difSum/counter;
    end
    
    %Moves Particles
    function [particles,xInitial,yInitial,theta] = moveParticles(p,x,y,xInt,yInt,clock)
        %Moves the particles in the direction given by their channel
        %INPUTS:
        %p - particle data (3 dimension array (xNodes by yNodes by 6))
        %x - x size of grid (Int)
        %y - y size of grid (Int)
        %xInt - initial x coordinate location of particles (3 dimension array (xNodes by yNodes by 6))
        %yInt - initial y coordinate location of particles (3 dimension array (xNodes by yNodes by 6))
        %clock - clock data (3 dimension array (xNodes by yNodes by 6))
        %OUTPUTS:
        %particles - new particle data where particles have moved (3 dimension array (xNodes by yNodes by 6))
        %xInitial - moved initial x coordinate data to follow the particles (3 dimension array (xNodes by yNodes by 6))
        %yInitial - moved initial y coordinate data to follow the particles (3 dimension array (xNodes by yNodes by 6))
        %theta - new clock data thats been moved (3 dimension array (xNodes by yNodes by 6))
    particles2=zeros(x,y,6);
    xInt2=zeros(x,y,6);
    yInt2 =zeros(x,y,6);
    clock2 = zeros(x,y,6);
        for i= 1:x
            for j = 1:y
                for k = 1:6
                    
                    if p(i,j,k) ==1
                       
                        %Move particle based on direciton while considering
                        %boundary conditions and update tracker if necessary
                        if k == 1
                            if(i+1 > x)
                                particles2(1,j,k)=1;
                                xInt2(1,j,k) = xInt(i,j,k);
                                yInt2(1,j,k)= yInt(i,j,k);
                                clock2(1,j,k)= clock(i,j,k);
                                
                            else
                                particles2(i+1,j,k)=1;
                                xInt2(i+1,j,k) = xInt(i,j,k);
                                yInt2(i+1,j,k)= yInt(i,j,k);
                                clock2(i+1,j,k)= clock(i,j,k);
                               
                            end
    
                        elseif k ==2
                            if(j+1>y)
                                particles2(i,1,k)= 1;
                                xInt2(i,1,k) = xInt(i,j,k);
                                yInt2(i,1,k)= yInt(i,j,k);
                                clock2(i,1,k)= clock(i,j,k);
                              
                            else
                                particles2(i,j+1,k)= 1;
                                xInt2(i,j+1,k) = xInt(i,j,k);
                                yInt2(i,j+1,k)= yInt(i,j,k);
                                clock2(i,j+1,k)= clock(i,j,k);
                            
                            end
    
                        elseif k ==3
                            if(i-1<1 &&j+1>y)
                                particles2(x,1,k)=1;
                                xInt2(x,1,k) = xInt(i,j,k);
                                yInt2(x,1,k)= yInt(i,j,k);
                                clock2(x,1,k)= clock(i,j,k);
                            elseif(i-1<1 && j+1 <=y)
                                particles2(x,j+1,k)=1;
                                xInt2(x,j+1,k) = xInt(i,j,k);
                                yInt2(x,j+1,k)= yInt(i,j,k);
                                clock2(x,j+1,k)= clock(i,j,k);
                            elseif(i-1>=1 && j+1 >y)
    
                                    particles2(i-1,1,k)=1;
                                    xInt2(i-1,1,k) = xInt(i,j,k);
                                    yInt2(i-1,1,k)= yInt(i,j,k);
                                    clock2(i-1,1,k)= clock(i,j,k);
                            else
                                particles2(i-1,j+1,k)= 1;
                                xInt2(i-1,j+1,k) = xInt(i,j,k);
                                yInt2(i-1,j+1,k)= yInt(i,j,k);
                                clock2(i-1,j+1,k)= clock(i,j,k);
                            end
    
                        elseif k ==4  
                            if(i-1<1)
                                particles2(x,j,k)=1;
                                xInt2(x,j,k) = xInt(i,j,k);
                                yInt2(x,j,k)= yInt(i,j,k);
                                clock2(x,j,k)= clock(i,j,k);
                            else
                                particles2(i-1,j,k)= 1;
                                xInt2(i-1,j,k) = xInt(i,j,k);
                                yInt2(i-1,j,k)= yInt(i,j,k);
                                clock2(i-1,j,k)= clock(i,j,k);
                            end
    
                        elseif k ==5
                            if(j-1<1)
                                    particles2(i,y,k)=1;
                                    xInt2(i,y,k) = xInt(i,j,k);
                                    yInt2(i,y,k)= yInt(i,j,k);
                                    clock2(i,y,k)= clock(i,j,k);

                            else
                                particles2(i,j-1,k)= 1;
                                xInt2(i,j-1,k) = xInt(i,j,k);
                                yInt2(i,j-1,k)= yInt(i,j,k);
                                clock2(i,j-1,k)= clock(i,j,k);
                            end
    
                        elseif k ==6
                            if(i+1>x &&j-1<1)
                                particles2(1,y,k)=1;
                                xInt2(1,y,k) = xInt(i,j,k);
                                yInt2(1,y,k)= yInt(i,j,k);
                                clock2(1,y,k)= clock(i,j,k);

                            elseif(i+1>x && j-1 >=1)
                                particles2(1,j-1,k)=1;
                                xInt2(1,j-1,k) = xInt(i,j,k);
                                yInt2(1,j-1,k)= yInt(i,j,k);
                                clock2(1,j-1,k)= clock(i,j,k);
                            elseif(i+1<=x && j-1 <1)
    
                                    particles2(i+1,y,k)=1;
                                    xInt2(i+1,y,k) = xInt(i,j,k);
                                    yInt2(i+1,y,k)= yInt(i,j,k);
                                    clock2(i+1,y,k)= clock(i,j,k);

                            else
                                particles2(i+1,j-1,k)= 1;
                                xInt2(i+1,j-1,k) = xInt(i,j,k);
                                yInt2(i+1,j-1,k)= yInt(i,j,k);
                                clock2(i+1,j-1,k)= clock(i,j,k);
                            end 
                        end
                    end
                end
            end
        end
        xInitial = xInt2;
        yInitial = yInt2;
        particles = particles2;
        theta = clock2;
    end
    
    %Interation step
    function [particles,theta]= ALGCAClocks2(p,x,y,t,J,J0,sticky)
        %Updates the configuration of the particles inside a hexagon
        %INPUTS:
        %p - particle data (3 dimension array (xNodes by yNodes by 6))
        %x - x gird size (Int)
        %y - y grid szie (Int)
        %t - clock data (3 dimension array (xNodes by yNodes by 6))
        %J - J parameter (Real)
        %J0 - J0 parameter (Real)
        %sticky - flag to know if using sticky rules (Boolean)
        %OUTPUTS:
        %particles - updated particle data (3 dimension array (xNodes by yNodes by 6))
        %theta - updated clock data (3 dimension array (xNodes by yNodes by 6))


    stencil = [1,0; 0,1;-1,1;-1,0;0,-1;1,-1;0,0;]; 
    carStencil=[1,0; cos(pi/3), sin(pi/3); -1 + cos(pi/3), sin(pi/3); -1,0; -cos(pi/3),-sin(pi/3); 1-cos(pi/3),-sin(pi/3)];
    %Calculates all permutations for 1,2,3,4,5,6
    perm1 = permutations(1);
    perm2 = permutations(2);
    perm3 = permutations(3);
    perm4 = permutations(4);
    perm5 = permutations(5);
    perm6 = permutations(6);
    
    p2 = zeros(x,y,6);
    t2 = zeros(x,y,6);
    
    for i =1:x
        for j = 1:y
            if sum(p(i,j,:))>0 %Test if hexagon has any particles
                %Calculates number of neighbors of a given particle taking
                %into consideration boundary conditions
                if(i-1<1 && j-1<1)
                    nNeighbors = sum(p(i+1,j,:)) + sum(p(i,j+1,:)) + sum(p(x,j+1,:)) +sum(p(x,j,:)) + sum(p(i,y,:)) + sum(p(i+1,y,:));
                elseif(i+1>x && j+1 >y)
                    nNeighbors = sum(p(1,j,:)) + sum(p(i,1,:)) + sum(p(i-1,1,:)) +sum(p(i-1,j,:)) + sum(p(i,j-1,:)) + sum(p(1,j-1,:));
                elseif(i+1>x && j-1<1)
                    nNeighbors = sum(p(1,j,:)) + sum(p(i,j+1,:)) + sum(p(i-1,j+1,:)) +sum(p(i-1,j,:)) + sum(p(i,y,:)) + sum(p(1,y,:));
                elseif(i-1<1 && j+1 > y)
                    nNeighbors = sum(p(i+1,j,:)) + sum(p(i,1,:)) + sum(p(x,1,:)) +sum(p(x,j,:)) + sum(p(i,j-1,:)) + sum(p(i+1,j-1,:));
                elseif(i+1> x)
                    nNeighbors = sum(p(1,j,:)) + sum(p(i,j+1,:)) + sum(p(i-1,j+1,:)) +sum(p(i-1,j,:)) +sum(p(i,j-1,:))+sum(p(1,j-1,:));
                elseif(i-1<1)
                    nNeighbors = sum(p(i+1,j,:)) + sum(p(i,j+1,:)) + sum(p(x,j+1,:)) +sum(p(x,j,:)) +sum(p(i,j-1,:))+sum(p(i+1,j-1,:));
                elseif(j+1>y)
                    nNeighbors = sum(p(i+1,j,:)) + sum(p(i,1,:)) + sum(p(i,1,:)) +sum(p(i-1,j,:)) +sum(p(i,j-1,:))+sum(p(i+1,j-1,:)); 
                elseif(j-1<1)
                    nNeighbors = sum(p(i+1,j,:)) + sum(p(i,j+1,:)) + sum(p(i-1,j+1,:)) +sum(p(i-1,j,:)) +sum(p(i,y,:))+sum(p(i+1,y,:));
                else
                    nNeighbors = sum(p(i+1,j,:)) + sum(p(i,j+1,:)) + sum(p(i-1,j+1,:)) +sum(p(i-1,j,:)) +sum(p(i,j-1,:))+sum(p(i+1,j-1,:));
                end
                A = zeros(6,nNeighbors,2); %A with hold what is basically the weighted gradient
                counter = 1;
                for k =1:6 %Iterate through each channel in a hexagon
                    if p(i,j,k) ~= 0
                       
                        % Calculating energy between target particle and its neighbors 
                        for in = 1:6 %Iterate through the neighboring hexagons 
                            neighborDirection = stencil(in,:); 
                            carNeighborDirection = carStencil(in,:);
                            for kn = 1:6 %Iterate through the channels of the neighboring hexagon
                                if(i+neighborDirection(1) > x && j+neighborDirection(2) <1) %Boundary Condition 
                                    if p(1,y,kn) ~= 0
                                        if(sticky)
                                            A(counter,6*(in-1)+kn,1)=neighborDirection(1)*(J0+J*cos(t(1,y,kn)));
                                            A(counter,6*(in-1)+kn,2)=neighborDirection(2)*(J0+J*cos(t(1,y,kn)));
                                        else
                                            A(counter,6*(in-1)+kn,1)=carNeighborDirection(1)*(J0+J*cos(t(1,y,kn)-t(i,j,k)));
                                            A(counter,6*(in-1)+kn,2)=carNeighborDirection(2)*(J0+J*cos(t(1,y,kn)-t(i,j,k)));
                                        end
    
                                    end
                                elseif(i+neighborDirection(1) < 1 && j+neighborDirection(2) >y)
                                    if p(x,1,kn) ~= 0
                                        if(sticky)
                                            A(counter,6*(in-1)+kn,1)=carNeighborDirection(1)*(J0+J*cos(t(x,1,kn)));
                                            A(counter,6*(in-1)+kn,2)=carNeighborDirection(2)*(J0+J*cos(t(x,1,kn)));
                                        else
                                            A(counter,6*(in-1)+kn,1)=carNeighborDirection(1)*(J0+J*cos(t(x,1,kn)-t(i,j,k)));
                                            A(counter,6*(in-1)+kn,2)=carNeighborDirection(2)*(J0+J*cos(t(x,1,kn)-t(i,j,k)));
                                        end
    
                                    end
                                elseif(i+neighborDirection(1) > x)
                                    if p(1,j+neighborDirection(2),kn) ~= 0
                                        if(sticky)
                                            A(counter,6*(in-1)+kn,1)=neighborDirection(1)*(J0+J*cos(t(1,j+neighborDirection(2),kn)));
                                            A(counter,6*(in-1)+kn,2)=neighborDirection(2)*(J0+J*cos(t(1,j+neighborDirection(2),kn)));
                                        else
                                            A(counter,6*(in-1)+kn,1)=carNeighborDirection(1)*(J0+J*cos(t(1,j+neighborDirection(2),kn)-t(i,j,k)));
                                            A(counter,6*(in-1)+kn,2)=carNeighborDirection(2)*(J0+J*cos(t(1,j+neighborDirection(2),kn)-t(i,j,k)));
                                        end
    
                                    end
                                elseif(i+neighborDirection(1) <1)
                                    if p(x,j+neighborDirection(2),kn) ~= 0
                                        if(sticky)
                                            A(counter,6*(in-1)+kn,1)=neighborDirection(1)*(J0+J*cos(t(x,j+neighborDirection(2),kn)));
                                            A(counter,6*(in-1)+kn,2)=neighborDirection(2)*(J0+J*cos(t(x,j+neighborDirection(2),kn)));
                                        else
                                            A(counter,6*(in-1)+kn,1)=carNeighborDirection(1)*(J0+J*cos(t(x,j+neighborDirection(2),kn)-t(i,j,k)));
                                            A(counter,6*(in-1)+kn,2)=carNeighborDirection(2)*(J0+J*cos(t(x,j+neighborDirection(2),kn)-t(i,j,k)));
                                        end
    
    
                                    end
                                elseif(j+neighborDirection(2) > y)
                                    if p(i+neighborDirection(1),1,kn) ~= 0
                                        if(sticky)
                                            A(counter,6*(in-1)+kn,1)=neighborDirection(1)*(J0+J*cos(t(i+neighborDirection(1),1,kn)));
                                            A(counter,6*(in-1)+kn,2)=neighborDirection(2)*(J0+J*cos(t(i+neighborDirection(1),1,kn)));
                                        else
                                            A(counter,6*(in-1)+kn,1)=carNeighborDirection(1)*(J0+J*cos(t(i+neighborDirection(1),1,kn)-t(i,j,k)));
                                            A(counter,6*(in-1)+kn,2)=carNeighborDirection(2)*(J0+J*cos(t(i+neighborDirection(1),1,kn)-t(i,j,k)));
                                        end
                                    end
                                    
                                elseif(j+neighborDirection(2) < 1)
    
                                        if p(i+neighborDirection(1),y,kn) ~= 0
                                            if(sticky)
                                                A(counter,6*(in-1)+kn,1)=neighborDirection(1)*(J0+J*cos(t(i+neighborDirection(1),y,kn)));
                                                A(counter,6*(in-1)+kn,2)=neighborDirection(2)*(J0+J*cos(t(i+neighborDirection(1),y,kn)));
                                            else
                                                A(counter,6*(in-1)+kn,1)=carNeighborDirection(1)*(J0+J*cos(t(i+neighborDirection(1),y,kn)-t(i,j,k)));
                                                A(counter,6*(in-1)+kn,2)=carNeighborDirection(2)*(J0+J*cos(t(i+neighborDirection(1),y,kn)-t(i,j,k)));
                                            end
    
                                        end
                                    
                                else
                                    if p(i+neighborDirection(1),j+neighborDirection(2),kn) ~= 0 %See if neighboring particle is there then calculate "gradient" with it
                                        if(sticky)
                                            A(counter,6*(in-1)+kn,1)=neighborDirection(1)*(J0+J*cos(t(i+neighborDirection(1),j+neighborDirection(2),kn)));
                                            A(counter,6*(in-1)+kn,2)=neighborDirection(2)*(J0+J*cos(t(i+neighborDirection(1),j+neighborDirection(2),kn)));
                                        else
                                            A(counter,6*(in-1)+kn,1)=carNeighborDirection(1)*(J0+J*cos(t(i+neighborDirection(1),j+neighborDirection(2),kn)-t(i,j,k)));
                                            A(counter,6*(in-1)+kn,2)=carNeighborDirection(2)*(J0+J*cos(t(i+neighborDirection(1),j+neighborDirection(2),kn)-t(i,j,k)));
                                        end
    
                                    end
                                end
    
                            end
                        end 
                        counter = counter +1;
                    end
                end

                %Calculates Gradient and flux 
                if(sum(p(i,j,:))==1) %Only looking at 1 permutations
                    transProbability = zeros(length(perm1),1);
                    
                    for iPerm = 1:length(perm1) %Calculates the transition probability 
                        gradient = [];
                        B=A;
    
                        flux = carStencil(iPerm,:);
                        B(1,:,1)= flux(1)*A(1,:,1);
                        B(1,:,2)=flux(2)*A(1,:,2);
                        gradient(1) = sum(sum(B(:,:,1)));
                        gradient(2) = sum(sum(B(:,:,2)));
                        
                        
                        transProbability(iPerm) = exp(gradient(1)+gradient(2));
                        
                    end
                    %Normalized probability
                    transProbability = transProbability/sum(transProbability);
                    %Randomly picks a permutation based on probabilities
                    condition = true;
                    ranNum = rand();
                    kk = 0;
                    sumProb = 0;
                    while(condition)
                        kk=kk+1;
                        sumProb = sumProb + transProbability(kk);
                        if(ranNum < sumProb)
                            condition = false;
                            index = kk;
                        end
                    end
                    %Assigns permutation
                    desiredPerm = perm1(index);
                    p2(i,j,desiredPerm) =1;
                    
                    %Updates maked particles
                    
                elseif(sum(p(i,j,:))==2) %Same thing as before but now for the 2 permutations
                    transProbability = zeros(length(perm2),1);
                    for iPerm = 1:length(perm2)
                        B=A;
                        gradient = [];

                        temp = perm2(iPerm,:);
                        channel1 = carStencil(temp(1),:);
                        channel2 = carStencil(temp(2),:);
                        B(1,:,1) = channel1(1)*A(1,:,1);
                        B(1,:,2) = channel1(2)*A(1,:,2);
                        B(2,:,1) = channel2(1)*A(2,:,1);
                        B(2,:,2) = channel2(2)*A(2,:,2);
                        gradient(1) = sum(sum(B(:,:,1)));
                        gradient(2) = sum(sum(B(:,:,2)));
                        
                        transProbability(iPerm) = exp(gradient(1)+gradient(2));
                        
                    end
                    transProbability = transProbability/sum(transProbability);
                    condition = true;
                    ranNum = rand();
                    kk = 0;
                    sumProb = 0;
                    while(condition)
                        kk=kk+1;
                        sumProb = sumProb + transProbability(kk);
                        if(ranNum < sumProb)
                            condition = false;
                            index = kk;
                        end
                    end
                    desiredPerm = perm2(index,:);
                    p2(i,j,desiredPerm(1)) =1;
                    p2(i,j,desiredPerm(2)) =1;
                    
                    
                elseif(sum(p(i,j,:))==3)%Same thing as before but now for the 3 permutations
                    transProbability = zeros(length(perm3),1);
                    for iPerm = 1:length(perm3)
                        B=A;
                        gradient = [];

                        temp = perm3(iPerm,:);
                        channel1 = carStencil(temp(1),:);
                        channel2 = carStencil(temp(2),:);
                        channel3 = carStencil(temp(3),:);
                        B(1,:,1) = channel1(1)*A(1,:,1);
                        B(1,:,2) = channel1(2)*A(1,:,2);
                        B(2,:,1) = channel2(1)*A(2,:,1);
                        B(2,:,2) = channel2(2)*A(2,:,2);
                        B(3,:,1) = channel3(1)*A(3,:,1);
                        B(3,:,2) = channel3(2)*A(3,:,2);
                        gradient(1) = sum(sum(B(:,:,1)));
                        gradient(2) = sum(sum(B(:,:,2)));
                        
                        transProbability(iPerm) = exp(gradient(1)+gradient(2));
                        
                    end
                    transProbability = transProbability/sum(transProbability);
                    condition = true;
                    ranNum = rand();
                    kk = 0;
                    sumProb = 0;
                    while(condition)
                        kk=kk+1;
                        sumProb = sumProb + transProbability(kk); 
                        if(ranNum < sumProb)
                            condition = false;
                            index = kk;
                        end
                    end
                    desiredPerm = perm3(index,:);
                    p2(i,j,desiredPerm(1)) =1;
                    p2(i,j,desiredPerm(2)) =1;
                    p2(i,j,desiredPerm(3)) =1;
    

                elseif(sum(p(i,j,:))==4) %Same thing as before but now for the 4 permutations
                    transProbability = zeros(length(perm4),1);
                    for iPerm = 1:length(perm4)
                        B=A;
                        gradient = [];

                        temp = perm4(iPerm,:);
                        channel1 = carStencil(temp(1),:);
                        channel2 = carStencil(temp(2),:);
                        channel3 = carStencil(temp(3),:);
                        channel4 = carStencil(temp(4),:);
                        B(1,:,1) = channel1(1)*A(1,:,1);
                        B(1,:,2) = channel1(2)*A(1,:,2);
                        B(2,:,1) = channel2(1)*A(2,:,1);
                        B(2,:,2) = channel2(2)*A(2,:,2);
                        B(3,:,1) = channel3(1)*A(3,:,1);
                        B(3,:,2) = channel3(2)*A(3,:,2);
                        B(4,:,1) = channel4(1)*A(4,:,1);
                        B(4,:,2) = channel4(2)*A(4,:,2);
                        gradient(1) = sum(sum(B(:,:,1)));
                        gradient(2) = sum(sum(B(:,:,2)));

                        transProbability(iPerm) = exp(gradient(1)+gradient(2));
                        
                    end
                    transProbability = transProbability/sum(transProbability);
                    condition = true;
                    ranNum = rand();
                    kk = 0;
                    sumProb = 0;
                    while(condition)
                        kk=kk+1;
                        sumProb = sumProb + transProbability(kk); 
                        if(ranNum < sumProb)
                            condition = false;
                            index = kk;
                        end
                    end
                    desiredPerm = perm4(index,:);
                    p2(i,j,desiredPerm(1)) =1;
                    p2(i,j,desiredPerm(2)) =1;
                    p2(i,j,desiredPerm(3)) =1;
                    p2(i,j,desiredPerm(4)) =1;
    
                elseif(sum(p(i,j,:))==5) %Same thing as before but now for the 5 permutations
                    transProbability = zeros(length(perm5),1);
                    for iPerm = 1:length(perm5)
                        B=A;
                        gradient = [];

                        temp = perm5(iPerm,:);
                        channel1 = carStencil(temp(1),:);
                        channel2 = carStencil(temp(2),:);
                        channel3 = carStencil(temp(3),:);
                        channel4 = carStencil(temp(4),:);
                        channel5 = carStencil(temp(5),:);
                        B(1,:,1) = channel1(1)*A(1,:,1);
                        B(1,:,2) = channel1(2)*A(1,:,2);
                        B(2,:,1) = channel2(1)*A(2,:,1);
                        B(2,:,2) = channel2(2)*A(2,:,2);
                        B(3,:,1) = channel3(1)*A(3,:,1);
                        B(3,:,2) = channel3(2)*A(3,:,2);
                        B(4,:,1) = channel4(1)*A(4,:,1);
                        B(4,:,2) = channel4(2)*A(4,:,2);
                        B(5,:,1) = channel5(1)*A(5,:,1);
                        B(5,:,2) = channel5(2)*A(5,:,2);
                        gradient(1) = sum(sum(B(:,:,1)));
                        gradient(2) = sum(sum(B(:,:,2)));

                        transProbability(iPerm) = exp(gradient(1)+gradient(2));
                        
                    end
                    transProbability = transProbability/sum(transProbability);
                    condition = true;
                    ranNum = rand();
                    kk = 0;
                    sumProb = 0;
                    
                    while(condition)
                        kk=kk+1;
                        sumProb = sumProb + transProbability(kk);
                        if(ranNum < sumProb)
                            condition = false;
                            index = kk;
                        end
                    end
                    desiredPerm = perm5(index,:);
                    p2(i,j,desiredPerm(1)) =1;
                    p2(i,j,desiredPerm(2)) =1;
                    p2(i,j,desiredPerm(3)) =1;
                    p2(i,j,desiredPerm(4)) =1;
                    p2(i,j,desiredPerm(5)) =1;
    
                elseif(sum(p(i,j,:))==6) %Same thing as before but now for the 6 permutations
                    hexR= abs(sum(exp(1i*t(i,j,:)))/6);
                    if(hexR <0.99)
                        transProbability = zeros(length(perm6),1);
                        for iPerm = 1:length(perm6)
                            B=A;
                            gradient = [];
    
                            temp = perm6(iPerm,:);
                            channel1 = carStencil(temp(1),:);
                            channel2 = carStencil(temp(2),:);
                            channel3 = carStencil(temp(3),:);
                            channel4 = carStencil(temp(4),:);
                            channel5 = carStencil(temp(5),:);
                            channel6 = carStencil(temp(6),:);
                            B(1,:,1) = channel1(1)*A(1,:,1);
                            B(1,:,2) = channel1(2)*A(1,:,2);
                            B(2,:,1) = channel2(1)*A(2,:,1);
                            B(2,:,2) = channel2(2)*A(2,:,2);
                            B(3,:,1) = channel3(1)*A(3,:,1);
                            B(3,:,2) = channel3(2)*A(3,:,2);
                            B(4,:,1) = channel4(1)*A(4,:,1);
                            B(4,:,2) = channel4(2)*A(4,:,2);
                            B(5,:,1) = channel5(1)*A(5,:,1);
                            B(5,:,2) = channel5(2)*A(5,:,2);
                            B(6,:,1) = channel6(1)*A(6,:,1);
                            B(6,:,2) = channel6(2)*A(6,:,2);
                            gradient(1) = sum(sum(B(:,:,1)));
                            gradient(2) = sum(sum(B(:,:,2)));

                            transProbability(iPerm) = exp(gradient(1)+gradient(2));
                            
                        end  
                        
                        
                        transProbability = transProbability/sum(transProbability);
                        condition = true;
                        ranNum = rand();
                        kk = 0;
                        sumProb = 0;
        
                        while(condition)
                            kk=kk+1;
                            sumProb = sumProb + transProbability(kk);
                            if(ranNum < sumProb)
                                condition = false;
                                index = kk;
                            end
                        end
                        desiredPerm = perm6(index,:);
                        p2(i,j,desiredPerm(1)) =1;
                        p2(i,j,desiredPerm(2)) =1;
                        p2(i,j,desiredPerm(3)) =1;
                        p2(i,j,desiredPerm(4)) =1;
                        p2(i,j,desiredPerm(5)) =1;
                        p2(i,j,desiredPerm(6)) =1;
 
                    else
                        p2(i,j,:)=p(i,j,:);
                        t2(i,j,:)=t(i,j,:);
        
                    end
                    
                end
                
                
                
            
            end
        end
    end
    particles = p2;
    theta = t2;
    
    end
    
    
    
    %Kuramoto Model
    function theta = kuramoto(p,x,y,theta,K,clockSpeed,delay)
        %Updates the clocks based on the kuramoto model
        %INPUTS:
        %p - particle data (3 dimension array (xNodes by yNodes by 6))
        %x - x grid size (Inr)
        %y - y grid size (Int)
        %theta - clocks data (3 dimension array (xNodes by yNodes by 6))
        %K - K Paramter (Real)
        %clockSpeed - the baseline speed of the clocks (Real)
        %delay - how delayed the clocks are in synchronizing (Real)
        %OUTPUTS: 
        %theta - the updated clock values (3 dimension array (xNodes by yNodes by 6))
        clock2 = theta;
        for i=1:x
            for j= 1:y
                for k = 1:6    
                    if(p(i,j,k)==1)
                        

                        sumKura = kuramotoSum(p,i,j,k,theta,x,y,delay); %Calculates the sum used in the kuramoto model
                        thetaChange = clockSpeed+K*sumKura; %calcualtes the change in the clock value
                       
    
                        clock2(i,j,k) = clock2(i,j,k) + thetaChange;
    
                    end
    
                end
            end
        end
        theta = clock2;
    end
    
    %Used in Kuramoto Model
    function sum = kuramotoSum(p,i,j,k,theta,xNodes,yNodes,delay)
        %Calculates the sum used in the Kuramoto model
        %INPUTS:
        %p - particle data (3 dimension array (xNodes by yNodes by 6))
        %i - the x index of the hexagon being examined (Int)
        %j - the y index of the hexagon being examined (Int)
        %theta - clock data (3 dimension array (xNodes by yNodes by 6))
        %xNodes - x grid size (Int)
        %yNodes - y grid Size (Int) 
        %delay - how delayed the clocks are in synchronizing (Real)
        stencil = [0,0;1,0; 0,1;-1,1;-1,0;0,-1;1,-1;];
        sum1 = 0;
        nNeighbours = 0;
        for iiKura = 2:length(stencil)
            temp = stencil(iiKura,:);
            

            for kk = 1:6 %Calculates the sum of sin(difference of clocks and neighbors) taking into account boundary conditions
                if(i+temp(1)<1 &&j+temp(2)>yNodes )
                    if(p(xNodes,1,kk) == 1)
                        nNeighbours = nNeighbours + 1;
                        sum1 = sum1 + sin(theta(xNodes,1,kk)-theta(i,j,k)-delay);
                    end
                elseif(i+temp(1)>xNodes && j+temp(2)<1)
                    if(p(1,yNodes,kk) == 1)
                        nNeighbours = nNeighbours + 1;
                        sum1 = sum1 + sin(theta(1,yNodes,kk)-theta(i,j,k)-delay);
                    end

                elseif(i+temp(1) >xNodes)
                    if(p(1,j+temp(2),kk) == 1)
                        nNeighbours = nNeighbours + 1;
                        sum1 = sum1 + sin(theta(1,j+temp(2),kk)-theta(i,j,k)-delay);
                    end
                elseif(j+temp(2)>yNodes)
                     if(p(i+temp(1),1,kk) == 1)
                        nNeighbours = nNeighbours + 1;
                        sum1 = sum1 + sin(theta(i+temp(1),1,kk)-theta(i,j,k)-delay);
                     end
                elseif(i+temp(1) <1)
                     if(p(xNodes,j+temp(2),kk) == 1)
                        nNeighbours = nNeighbours + 1;
                        sum1 = sum1 + sin(theta(xNodes,j+temp(2),kk)-theta(i,j,k)-delay);
                     end
                elseif(j+temp(2)<1)
                     if(p(i+temp(1),yNodes,kk) == 1)
                        nNeighbours = nNeighbours + 1;
                        sum1 = sum1 + sin(theta(i+temp(1),yNodes,kk)-theta(i,j,k)-delay);
                     end  
                else
                    if(p(i+temp(1),j+temp(2),kk) == 1)
                        nNeighbours = nNeighbours + 1;
                        sum1 = sum1 + sin(theta(i+temp(1),j+temp(2),kk)-theta(i,j,k)-delay);
                    end
                end
            end

    
        end
        if(nNeighbours ~= 0)
            sum =sum1/nNeighbours;
        else
            sum = 0;
        end
    end
    
    % R value test
    function rMag =radius(p,x,y,theta)
        %Calculates the r value of the simulation
        %INPUTS: 
        %p - particle data (3 dimension array (xNodes by yNodes by 6))
        %x - x grid size (Inr)
        %y - y grid size (Int)
        %theta - clocks data (3 dimension array (xNodes by yNodes by 6))
        %OUTPUTS: 
        %rMag - R Value
        sumR = 0;
        counter = 0;
        for i = 1:x
            for j = 1:y
                for k = 1:6
                    if p(i,j,k)==1
                        sumR = sumR +exp(1i*theta(i,j,k));
                        counter = counter +1; 
                    end
                end
            end
        end
        r = sumR/counter;
        rMag = abs(r);
    
    end
    
    %Generates unique permutations of V
    function P = uniqueperms(V)
        %INPUTS: 
        %V - array of numbers (array)
        %OUTPUTS:
        %P - array of all possible permutations (array)
        NPerms = numel(V);
        P = reshape(V,1,NPerms);
        if NPerms < 2 
            return; 
        end
         
            [W, ~, IX] = unique(V);
        if numel(W) == 1
            return;
        end
        IX = sort(reshape(IX, 1, NPerms));
        KPerms = histcounts(IX, 1:IX(end)+1);
        FPerms= cumprod([1 1 2:NPerms]);
        MPerms = FPerms(NPerms+1) / prod(FPerms(KPerms+1));
        P = repmat(IX,MPerms,1);
        Q = P(1,:);
        for k = 2:MPerms 
            i = find(Q(2:end) > Q(1:end-1),1,'last');
            j = find(Q(i) < Q,1,'last');
            Q([i j]) = Q([j i]) ;
            Q(i+1:NPerms) = Q(NPerms:-1:i+1);
            P(k,:) = Q;
        end
        P = W(P);
    end
    
    %Generates all possible permuations of size nPerms using 1,2,3,4,5,6
    function perms = permutations(nPerms)
        %INPUT:
        %nPerms - size of the permutation
        %OUTPUT: 
        %perms - all possible permutations of size nPerms that consist only of 1,2,3,4,5,6
        if(nPerms == 1)
            perms = [1,2,3,4,5,6];
        elseif(nPerms ==2)
            perms=zeros(30,2);
            permCount=1;
            for i = 1:6
                for j = 1:6
                    if (i~= j)
                        perms(permCount,1) =i;
                        perms(permCount,2) =j;
                        permCount = permCount +1;
                    end
                end
            end
        elseif(nPerms==3)
            perms=zeros(120,3);
            permCount=1;
            for i = 1:6
                for j =1:6
                    for k = 1:6
                        if(i~=j)
                            if(i~=k)
                                if(j~=k)
                                   perms(permCount,1) =i;
                                   perms(permCount,2) =j;
                                   perms(permCount,3) =k;
                                   permCount = permCount +1; 
                                end
                            end
                        end
                    end
                end
            end 
        elseif(nPerms==4)
            perms=zeros(360,4);
            permCount=1;
            for i = 1:6
                for j =1:6
                    for k = 1:6
                        for l =1:6
                            if(i~=j)
                                if(i~=k)
                                    if(i~=l)
                                        
                                        if(j~=k)
                                            if(j~=l)
                                                if(k~=l)
                                                   perms(permCount,1) =i;
                                                   perms(permCount,2) =j;
                                                   perms(permCount,3) =k;
                                                   perms(permCount,4) =l;
                                                   permCount = permCount +1; 
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end    
        
        elseif(nPerms==5)
            perms=zeros(720,5);
            permCount=1;
            for i = 1:6
                for j =1:6
                    for k = 1:6
                        for l =1:6
                            for m=1:6
                                if(i~=j)
                                    if(i~=k)
                                        if(i~=l)
                                            if(i~=m)
        
                                                if(j~=k)
                                                    if(j~=l)
                                                        if(j~= m)
                                                            
                                                            if(k~=l)
                                                                if(k~=m)
                                                                    if(l~=m)
                                                                       perms(permCount,1) =i;
                                                                       perms(permCount,2) =j;
                                                                       perms(permCount,3) =k;
                                                                       perms(permCount,4) =l;
                                                                       perms(permCount,5) =m;
                                                                       permCount = permCount +1; 
                                                                    end
                                                                
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        elseif(nPerms==6)
            perms = uniqueperms([1,2,3,4,5,6]);  
        end
    end
    
    
   
    % function nCount = neighborhoodCount(size)
    %     if size ==0
    %         nCount = 1;
    %     else
    %         nCount = size*6 + neighborhoodCount(size-1);
    %     end
    % end
    
    function avgLocalR = localRvalue(theta,x,y)
        %Calculates the r value on a local level (neighborhood of size 1)
        %INPUTS:
        %theta- clock data
        %x- x grid size
        %y - y grid size
        stencil = [0,0;1,0; 0,1;-1,1;-1,0;0,-1;1,-1;];
        counter = 0;

        for i = 1:x
            for j = 1:y
                if sum(theta(i,j,:))>0
                    for k = 1:6
                       if theta(i,j,k)>0
                           
                           counter = counter+1;
                           nNeighbours =0;
                           sum1 = 0;
                           for iiLocalR = 1:length(stencil)
                                temp = stencil(iiLocalR,:);
                                for kk = 1:6
                                    if(i+temp(1)<1 &&j+temp(2)>y )
                                        if(theta(x,1,kk)> 0)
                                            nNeighbours = nNeighbours + 1;
                                            sum1 = sum1 + exp(1i*theta(x,1,kk)); 
                                        end
                                    elseif(i+temp(1)>x && temp(2)<1)
                                        if(theta(1,y,kk)> 0)
                                            nNeighbours = nNeighbours + 1;
                                            sum1 = sum1 + exp(1i*theta(1,y,kk)); 
                                        end
                                    elseif(i+temp(1) >x)
                                        if(theta(1,j+temp(2),kk)>0)
                                            nNeighbours = nNeighbours + 1;
                                            sum1 = sum1 + exp(1i*theta(1,j+temp(2),kk));
                                        end
                                    elseif(j+temp(2)>y)
                                         if(theta(i+temp(1),1,kk)> 0)
                                            nNeighbours = nNeighbours + 1;
                                            sum1 = sum1 + exp(1i*theta(i+temp(1),1,kk)); 
                                         end
                                    elseif(i+temp(1) <1)
                                         if(theta(x,j+temp(2),kk) >0)
                                            nNeighbours = nNeighbours + 1;
                                            sum1 = sum1 + exp(1i*theta(x,j+temp(2),kk));
                                         end
                                    elseif(j+temp(2)<1)
                                         if(theta(i+temp(1),y,kk) >0)
                                            nNeighbours = nNeighbours + 1;
                                            sum1 = sum1 + exp(1i*theta(i+temp(1),y,kk)); 
                                         end     
                                    else
                                        if(theta(i+temp(1),j+temp(2),kk) >0)
                                            nNeighbours = nNeighbours + 1;
                                            sum1 = sum1 + exp(1i*theta(i+temp(1),j+temp(2),kk)); 
                                        end
                                    end
                                end
                           end
                           localRFunc(counter) = abs(sum1/nNeighbours);
                       end
                    end
                end
            end
        end
        avgLocalR = mean(localRFunc);
    
    end
    function [aggregates] = cleanAggregates(x,y,particles)
        %Cleans the aggregates by removing particles that are in hexagons
        %that have a total of 2 or less
        %INPUTS: 
        %x - x grid size
        %y - y grid size
        %particles - particle data
        %OUTPUTS:
        %aggregates - cleaned up particle data
        aggregates = zeros(x,y);
            for i = 1:x
                for j = 1:y
                    if sum(particles(i,j,:))>2
  
                        for k= 1:6
                            if particles(i,j,k)==1

                                aggregates(i,j,k)=1;
                            end
                        end
                    end
                end
            end
    end
    function plotParticlesSimple(p,x,y,dotSize)
        %Same has PlotParticles2 but dots are black
        %INPUTS: 
        %p - particle data
        %x - x grid size
        %y - y grid size
        %dotSize - size of dots being plotted
    stencil = [(1/3),0; 0,(1/3);-(1/3),(1/3);-(1/3),0;0,-(1/3);(1/3),-(1/3);0,0;]; 
    Points = zeros(sum(sum(sum(p))),2);
    pointCounter = 1;
        for i =1:x
            for j = 1:y
                for k = 1:6
                    if p(i,j,k)==1
                        
                        %Calculating where each particle needs to be placed and
                        temp = stencil(k,:);
                        if(i>y-floor(j/2)) %i<ceil((y-j)/2))%Changing the shape of the layout to a rectangle for display purposes %%%i>y-floor(j/2)
                            Points(pointCounter,1)= (i-x)+temp(1)+ (j+temp(2))*cos(pi/3);
                            Points(pointCounter,2)= (j+temp(2))*sin(pi/3);
    
                        else
                            
                            Points(pointCounter,1)= i+temp(1)+ (j+temp(2))*cos(pi/3);
                            Points(pointCounter,2)= (j+temp(2))*sin(pi/3);
    
                            
                        end
                        pointCounter = pointCounter + 1;
        
                    end
                end
            end
        end
        
        figure(f2)
        clf(2)
        axis off
        scatter(Points(:,1),Points(:,2),dotSize,'black','filled');
        f = gcf;
        f.Position = [50 50 1000 725];
        set(gcf,'color','w');
        
     
        xlim([0 x+1])
        ylim([0 y])
        set(gca,'LooseInset',get(gca,'TightInset'));
        axis off
    end
end
function [aggCount] = countAggregates(imageX)
    %Counts the number of aggregates of the image
    %INPUTS:
    %imageX - image of aggregates going to be counted
    %OUTPUTS:
    %aggCount - the number of aggregates in the image
    imageX = bwareaopen(imageX,40); %fills in holes of the aggregates

    %the next 50 or so lines are used to erase white space along the edge
    %of the image. It does this by deleting any row/col that is entirely
    %zeros 
    findZeros = true;
    count = 1;
    while(findZeros)
        if(isempty(find(~imageX(count,:),1))==0)
            if(count~=1)
                imageX(1:count-1,:)=[];
                findZeros = false;
            else
                findZeros=false;
            end
        end
        count = count+1;
    end
    findZeros = true;
    count = 1;
    while(findZeros)
        if(isempty(find(~imageX(:,count),1))==0)
            if(count~=1)
                imageX(:,1:count-1)=[];
                findZeros = false;
            else
                findZeros=false;
            end
        end
        count = count +1;
    end
    findZeros = true;
    [M N] = size(imageX);
    M2=M;
    N2=N;
    while(findZeros)
        if(isempty(find(~imageX(M2,:),1))==0)
            if(M2~=M)
                imageX(M2+1:M,:)=[];
                findZeros = false;
            else
                findZeros=false;
            end
        end
        M2 = M2-1;
    end
    findZeros = true;
    while(findZeros)
        if(isempty(find(~imageX(:,N2),1))==0)
            if(N2~=N)
                imageX(:,N2+1:N)=[];
                findZeros = false;
            else
                findZeros=false;
            end
        end
        N2 = N2-1;
    end
    [M N] = size(imageX);

    %Create a big image by combining 91 of the original sized image
    bigX = [ones(M,floor(2.5*N)), imageX,imageX,imageX,imageX,imageX,imageX,ones(M,ceil(2.5*N));...
        ones(M,2*N),imageX,imageX,imageX,imageX,imageX,imageX,imageX,ones(M,2*N);...
        ones(M,floor(1.5*N)),imageX,imageX,imageX,imageX,imageX,imageX,imageX,imageX,ones(M,ceil(1.5*N));...
        ones(M,N),imageX,imageX,imageX,imageX,imageX,imageX,imageX,imageX,imageX,ones(M,N);...
        ones(M,floor(N/2)),imageX,imageX,imageX,imageX,imageX,imageX,imageX,imageX,imageX,imageX,ones(M,ceil(N/2));...
        imageX,imageX,imageX,imageX,imageX,imageX,imageX,imageX,imageX,imageX,imageX;...
        ones(M,floor(N/2)),imageX,imageX,imageX,imageX,imageX,imageX,imageX,imageX,imageX,imageX,ones(M,ceil(N/2));...
        ones(M,N),imageX,imageX,imageX,imageX,imageX,imageX,imageX,imageX,imageX,ones(M,N);...
        ones(M,floor(1.5*N)),imageX,imageX,imageX,imageX,imageX,imageX,imageX,imageX,ones(M,ceil(1.5*N));...
        ones(M,2*N),imageX,imageX,imageX,imageX,imageX,imageX,imageX,ones(M,2*N);...
        ones(M,floor(2.5*N)), imageX,imageX,imageX,imageX,imageX,imageX,ones(M,ceil(2.5*N))];
    
    %Smooth up the edges
    windowSize=12; 
    kernel=ones(windowSize)/windowSize^2;
    result=conv2(single(bigX),kernel,'same');
    result=result>0.5;
    bigX(~result)=0;

    CC = bwconncomp(~bigX);
    aggCount =floor(CC.NumObjects/91);
end
function [avgRadius,avgCirc,stdArea] = radiusOfG(X)
    %Calculates radius of gyration, circularity and std area given an image
    %INPUTS:
    %X - image
    %avgRadius - the average radius of gyration
    %avgCirc - average circularity
    %stdArea - the standard deviation of the area of the aggregates
    
    X2= bwareaopen(X,100); %fills in holes in the image

    %Eliminates white space along the edge
    findZeros = true;
    count = 1;
    while(findZeros)
        if(isempty(find(~X2(count,:),1))==0)
            if(count~=1)
                X2(1:count-1,:)=[];
                findZeros = false;
            else
                findZeros=false;
            end
        end
        count = count+1;
    end
    findZeros = true;
    count = 1;
    while(findZeros)
        if(isempty(find(~X2(:,count),1))==0)
            if(count~=1)
                X2(:,1:count-1)=[];
                findZeros = false;
            else
                findZeros=false;
            end
        end
        count = count +1;
    end
    findZeros = true;
    [M N] = size(X2);
    M2=M;
    N2=N;
    while(findZeros)
        if(isempty(find(~X2(M2,:),1))==0)
            if(M2~=M)
                X2(M2+1:M,:)=[];
                findZeros = false;
            else
                findZeros=false;
            end
        end
        M2 = M2-1;
    end
    findZeros = true;
    while(findZeros)
        if(isempty(find(~X2(:,N2),1))==0)
            if(N2~=N)
                X2(:,N2+1:N)=[];
                findZeros = false;
            else
                findZeros=false;
            end
        end
        N2 = N2-1;
    end
    [M N] = size(X2);
    
    %Create big image of 91 original images
    bigX = [ones(M,floor(2.5*N)), X2,X2,X2,X2,X2,X2,ones(M,ceil(2.5*N));...
        ones(M,2*N),X2,X2,X2,X2,X2,X2,X2,ones(M,2*N);...
        ones(M,floor(1.5*N)),X2,X2,X2,X2,X2,X2,X2,X2,ones(M,ceil(1.5*N));...
        ones(M,N),X2,X2,X2,X2,X2,X2,X2,X2,X2,ones(M,N);...
        ones(M,floor(N/2)),X2,X2,X2,X2,X2,X2,X2,X2,X2,X2,ones(M,ceil(N/2));...
        X2,X2,X2,X2,X2,X2,X2,X2,X2,X2,X2;...
        ones(M,floor(N/2)),X2,X2,X2,X2,X2,X2,X2,X2,X2,X2,ones(M,ceil(N/2));...
        ones(M,N),X2,X2,X2,X2,X2,X2,X2,X2,X2,ones(M,N);...
        ones(M,floor(1.5*N)),X2,X2,X2,X2,X2,X2,X2,X2,ones(M,ceil(1.5*N));...
        ones(M,2*N),X2,X2,X2,X2,X2,X2,X2,ones(M,2*N);...
        ones(M,floor(2.5*N)), X2,X2,X2,X2,X2,X2,ones(M,ceil(2.5*N))];
    
    %smooth edges
    kernel = ones(25)/25^2;
    Y = convolve2(double(bigX),kernel,'wrap');
    bigX = Y>0.5;

    [M, N] = size(bigX);
    CC=bwconncomp(~bigX);
    S = regionprops(CC,'Centroid','Circularity','Area');
    centroids = cat(1,S.Centroid);
    circularity = cat(1,S.Circularity);
    area=cat(1,S.Area);
    numOb = CC.NumObjects;
    radius = zeros(numOb,1);

    for i = 1:numOb %Calculates the radius of gyration of each aggregate
        x = CC.PixelIdxList{1,i};
        [row, col] = ind2sub([M, N], x);
        xCent = centroids(i,1);
        yCent = centroids(i,2);
        norms = zeros(length(row),1);

        for j=1:length(row)
            norms(j) = (norm([col(j), row(j)] - [xCent,yCent]))^2;  
        end
        radius(i) = (sqrt(sum(norms)/length(row)));




    end

    avgRadius = mean(radius);
    avgCirc = mean(circularity);
    stdArea = std(area,1);
    
end

