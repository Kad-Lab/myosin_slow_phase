clf;
%% Iterate parameters within the model
clearvars;
loops = 1;
chunkID = 1; force = 0;
% Set output folder to wherever this script is being run from
outDir = pwd;   % captures the folder you navigated to before hitting Run
%Uncomment to see model
%     figure
%     I = imread("Model.jpg"); 
%     imshow(I);
    
%top traverse of scheme
kf1 = 0.01; kr1 = 0.51;         %MB-MB=MD-MB
kf3 = 0.01; kr3 = 0.35;         %MD-MB=MD-MD
%bottom traverse of scheme
kf2 = 0; kr2 = 0;       %MB-MB=MB-MD
kf4 = 0; kr4 = 0;       %MB-MD=MD-MD

loops = 1; 
loopmax = 10000; %number of MD2_ATP values - check line 294 as well

%For running with iteration of all parameters
for kf3 = 0.01:0.05:0.51
    for kr3 = 0.01:0.05:0.51 %to nest a loop of equ 3 into a variation of equ 1
        for kf1 = 0.01:0.05:0.51 %to nest a loop of equ 3 into a variation of equ 1
            for kr1 = 0.01:0.05:0.51
                [MD1_ATP, MD2_ATP] = Main_Cycle(kf1, kr1, kf2, kr2, kf3, kr3, kf4, kr4, loopmax);
                MDx_ATP = [MD1_ATP; MD2_ATP]; %puts MD1_ATP and MD2_ATP together
                MDx_ATP = MDx_ATP(MDx_ATP~=0); %Clears out zeros
                MD1_ATP_Stor(1:size(MDx_ATP,1)+8,loops) = [kf1; kr1; kf2; kr2; kf3; kr3; kf4; kr4; MDx_ATP];
                %MD2_ATP_Stor(1:size(MD2_ATP,1)+8,loops) = [kf1; kr1; kf2; kr2; kf3; kr3; kf4; kr4; MD2_ATP];
                loops = loops +1; 
                % size check + flush happens here
                [MD1_ATP_Stor, loops, chunkID, force] = flushIfTooBig(MD1_ATP_Stor, loops, chunkID, force, outDir);
            end
        end
    end
end

% name = ['MD1_Stor_kr2_' num2str(kf3) '_.mat'];
% save(name,'MD1_ATP_Stor');
% name = ['MD2_Stor_kr2_' num2str(kf3) '_.mat'];
% save(name,'MD2_ATP_Stor');

vName = sprintf('MD1_%04d', chunkID);

% Save (create file with new name)
save(fullfile(outDir, vName), 'MD1_ATP_Stor', '-v7.3');
% Measure MD1 size in RAM
w = whos('MD1_ATP_Stor');
bytes = w.bytes;

fprintf('Flushed %s (%.2f GB) to %s\n', vName, bytes/1e9, vName);

%---------------------------------------

% %If you want to run this without iteration
% for i = 1:loops %to run this multiple times
%     [MD1_ATP, MD2_ATP, myostate] = Main_Cycle(kf1, kr1, kf2, kr2, kf3, kr3, kf4, kr4, loopmax);
%     %[MD1_ATP,MD2_ATP] = fit_data(MD1_ATP,MD2_ATP);
%     MD_stor(i,:) = [MD1_ATP(2) MD2_ATP(2)];
% end
% Save (create file with new name)
% save(fullfile(outDir, vName), 'MD1_ATP_Stor', '-v7.3');
% Measure MD1 size in RAM
% w = whos('MD1_ATP_Stor');
% bytes = w.bytes;
% [file, location] = uiputfile('*.xlsx', 'Save As');
% if file ~= 0
%     fullPath = fullfile(location, file);
%     writematrix(MD_stor, fullPath);
% end

%Show curvefits
% [MD1_ATP,MD2_ATP,fitted_curve] = fit_data(MD1_ATP,MD2_ATP); %go fit the data and in the next version save the fitted_curve somehow with the rates
% disp('MD1/MD2 ratio:');
% size(MD1_ATP)/size(MD2_ATP)


%% Flush variables function
function [MD1_ATP_Stor, loops, chunkID,force] = flushIfTooBig(MD1_ATP_Stor, loops, chunkID, force, outDir)
% flushIfTooBig
% Flushes a single large variable (MD1_ATP_Stor) to disk in chunks when it exceeds a threshold.
%
% Variables created on disk: MD1_0001, MD1_0002, ...
%
% Inputs:
%   MD1     - your in-memory buffer
%   loops   - your column counter (will be reset to 1 after flush)
%   chunkID - chunk counter you keep in the main loop
%   loops is used to increment the column for MD1_ATP_Stor 
% Outputs:
%   MD1     - emptied after flush
%   loops   - reset to 1 after flush
%   chunkID - incremented after flush



    % ---- SETTINGS ----

    maxBytes = 5e9;  % 5 GB threshold (change to required)
    % ------------------

    if isempty(MD1_ATP_Stor)
        return
    end

    % Measure MD1 size in RAM
    w = whos('MD1_ATP_Stor');
    bytes = w.bytes;

    if bytes < maxBytes
        return
    end

    % Create chunk variable name
    vName = sprintf('MD1_%04d', chunkID);

    % Save (create file with new name)
    save(fullfile(outDir, vName), 'MD1_ATP_Stor', '-v7.3');

    fprintf('Flushed %s (%.2f GB) to %s\n', vName, bytes/1e9, vName);

    % Reset in-memory buffer + counters
    MD1_ATP_Stor = [];
    loops = 1;
    chunkID = chunkID + 1;

end



  %% Fit data function
function [MD1_ATP,MD2_ATP,fitted_curve] = fit_data(MD1_ATP,MD2_ATP)
%% Fitting routine as a sub-section in the function use CTRL+Enter to run the fit
clearvars h;
% MD1_stor=nonzeros(MD1_stor);
% MD2_stor=nonzeros(MD2_stor);
% total_ATP = cat(1,MD1_stor,MD2_stor);
% histogram(total_ATP);
% histogram(cat(1,MD1_ATP,MD2_ATP));

% [N, edges] = histcounts(MD1_stor);%,'normalization','probability');
% centers = (edges(1:end-1) + edges(2:end))./2; %histcounts gives the bin edges, but we want to plot the bin centers
% plot(centers,N,'ko');
% set(gca, 'YScale', 'log')

g = cat(1,MD1_ATP,MD2_ATP); %Put both vectors together end-to-end
%  g = MD2_ATP;

% [N, edges] = histcounts(g); % Get the histogram values
% N = cat(2,N,0); %address shorter N than edges
% h = [edges' log(N)']; %create a new vector with two columns of the histogram
% h = h( ~any(isnan(h)|isinf(h)|h==0,2 ),: );  

%Create a cumulative frequency histogram
time_range = max(g); %Can replace with max(g) or 500
no_bins = 500; %number of bins
h(:,1) = linspace(0, time_range, no_bins); 
for i = 1:size(h(:,1),1)
    h(i,2) = sum(g>h(i,1) | g==h(i,1)); %Put values into col2
    if i>2 && h(i-1,2) == h(i,2) %Remove duplicates
        h(i,2) = inf;
    end
end
h(:,2) = log(h(:,2));
f = h(1,:);
h = h( ~any(isnan(h)|isinf(h)|h==0,2 ),: ); 
h = cat(1,f,h);

plot(h(:,1),(h(:,2)),'ko');


%single exponentials
% fitfun = fittype(@(k1,a1,x) log(exp(-k1*x)*a1));
% [fitted_curve, gof] = fit(h(:,1),h(:,2),fitfun,'StartPoint', [0.05 numel(g)/100]);
% coeffvals = coeffvalues(fitted_curve); 
% 
% fitted_curve %Show fit on command line
% hold on;
% plot(h(:,1),fitted_curve(h(:,1)),'DisplayName','single'); %Plot fit for single exp

%Two exponentials
fitfun = fittype(@(k1,a1,k2,a2,x) log(exp(-k1*x)*a1+exp(-k2*x)*a2));
[fitted_curve, gof] = fit(h(:,1),h(:,2),fitfun,'StartPoint', [0.05 numel(g) 0.01 numel(g)/100]);
coeffvals = coeffvalues(fitted_curve); 

fitted_curve %Show fit on command line
hold on;
plot(h(:,1),fitted_curve(h(:,1)),'DisplayName','double'); %Plot fit for double exp
% disp('a1/a2 ratio:');
% coeffvals(1,2)/coeffvals(1,4) %Display ratios of amplitudes
disp('%fast:');
100*coeffvals(1,2)/(coeffvals(1,2)+coeffvals(1,4)) %Display percentage of fast amplitude

%Three exponentials
% fitfun = fittype(@(k1,a1,k2,a2,k3, a3,x) log(exp(-k1*x)*a1+exp(-k2*x)*a2+exp(-k3*x)*a3));
% [fitted_curve, gof] = fit(h(:,1),h(:,2),fitfun,'StartPoint', [0.05 numel(g)/2 0.0017 numel(g)/10 0.0001 numel(g)/10]);
% coeffvals = coeffvalues(fitted_curve);
% 
% 
% fitted_curve %Show fit on command line
% plot(h(:,1),fitted_curve(h(:,1)),'DisplayName','triple'); %Plot fit for triple exp

% axis([0, 1000 0.98*min(h(:,2)) 1.02*max(h(:,2))]); %use to change axis limits
legend ;
hold off;


end


 %% Function to determine if ATP is released after coming from an active state
function [MD1,MD2,MD1_ATP,MD2_ATP] = ATP_rel(kATP,MD1,MD2,percent_ATP,MD1_ATP,MD2_ATP)
%MDx includes 1 = fluorescent state, 2 = ATP_life, 3 = MDx_life, 4 = Total ATP lifetime
flag1 = 0; flag2 = 0; MD = 1; %Currently MDx(1) is MD and fixed to 1, so is a spare register.
rebound1 = false; rebound2 = false;   % true if the current ATP was bound inside this call

    while or(flag1 == 0, flag2 == 0)
        %first head
        
        if MD1(2)<MD1(3) %is ATP lifetime less than the total head lifetime?
            leftover_ATP_life = MD1(3)-MD1(2); %This removes the ATP lifetime from the head lifetime to give an accurate remaining lifetime
%             [MD,ATP_life] = ATP_lifetime(kATP,percent_ATP); %determine new ATP lifetime and if fluorescent
                if rebound1
                    logged1 = MD1(4);                        % full age of rebound ATP
                else
                    logged1 = MD1(4)-leftover_ATP_life;      % ATP present at state entry
                end

            ATP_life = -log(rand(1))/(kATP); %determine new ATP lifetime
            MD1_ATP = cat(1,MD1_ATP,logged1); 
            MD1 = [MD,ATP_life,leftover_ATP_life,leftover_ATP_life]; %rebuild the head 1 state
            if MD1(2)<MD1(3)
                flag1 = 0; rebound1 = true;
                MD1(4) = ATP_life;
            else
                flag1 = 1; rebound1 = false;
            end
        else
                flag1 = 1;
        end
        
        
        % second head
        
        if MD2(2)<MD2(3) %is ATP lifetime less than the total head lifetime?
            leftover_ATP_life = MD2(3)-MD2(2); %This removes the ATP lifetime from the head lifetime to give an accurate remaining lifetime
%             [MD,ATP_life] = ATP_lifetime(kATP,percent_ATP); %determine new ATP lifetime and if fluorescent
                if rebound2
                    logged2 = MD2(4);
                else
                    logged2 = MD2(4)-leftover_ATP_life;
                end
            MD2_ATP = cat(1,MD2_ATP,logged2);


            ATP_life = -log(rand(1))/(kATP); %determine new ATP lifetime
            
            MD2 = [MD,ATP_life,leftover_ATP_life,leftover_ATP_life];%MD2 = [MD,ATP_life,0,0]; %rebuild the head 2 state
            if MD2(2)<MD2(3)
                flag2 = 0; rebound2 = true;
                MD2(4) = ATP_life;
            else
                flag2 = 1; rebound2 = false;
            end
        else
                flag2 = 1;
        end
    
    end

end
%% Main cycle function
function [MD1_ATP, MD2_ATP, myostate] = Main_Cycle(kf1, kr1, kf2, kr2, kf3, kr3, kf4, kr4, loopmax)
            % 'Monte Carlo Simulation of the IHM two step model (MAG) model' ...
%     

% Setup overall loops


myosin = 1; %Molecular state 1=MBMB, 2=MDMB, 3=MBMD, 4=MDMD
myolife = 0; life = 0;
myostate = [0,0,0,0]; %Count for number of visitations
MD = 0; MD1(4) = 0; MD2(4) = 0; %Clear out the heads to begin with (four columns first is ATP present? Second ATP_life, third head lifetime, fourth total ATP lifetime)
MD1_life = 0; MD2_life = 0; %For calculation of lifetimes after the simulation
MD1_ATP = 0; MD2_ATP = 0; %These arrays carry the lifetimes of each ATP from MDx(4)

kATP = 0.05; %0.05/s ATP turnover time. Set at 0.5/s for Pi release being limiting (see notebook 8/1/26)
percent_ATP = 1; %percent amount of fluorescent ATP 


while size(MD1_ATP,1) < loopmax || size(MD2_ATP,1) < loopmax %Cycle until enough data is collected  ,size(MD2_ATP,1)) 
% Then check to see what state myosin is in and give it a
% lifetime and then decide where it goes.
    if myosin == 1
        MD1(3) = 0; %This resets the lifetimes for head 1 because it is down in state 1
        MD2(3) = 0; %This resets the lifetimes for head 2 because it is down in state 1
        life = -log(rand(1))/(kf1+kf2); %give the state a lifetime
        MD1(4) = MD1(4) + life; %keep track of total time before ATP is released
        MD2(4) = MD2(4) + life; %keep track of total time before ATP is released
        myolife = myolife + life;
        myostate(1) = myostate(1) + 1;
        if rand(1)< (kf1/(kf1+kf2))
            myosin = 2;
            life = -log(rand(1))/(kr1+kf3); %give the state a lifetime
            MD1(2) = -log(rand(1))/(kATP); %From state 1 so refresh ATP lifetime
            MD1(3) = MD1(3) + life; %add state2 lifetime to previous lifetime history
            MD1(4) = MD1(4) + life; %keep track of total time before ATP is released
            MD2(4) = MD2(4) + life; %keep track of total time before ATP is released
            myolife = myolife + life; %myolife is the total clock progression
            [MD1,MD2,MD1_ATP,MD2_ATP] = ATP_rel(kATP,MD1,MD2,percent_ATP,MD1_ATP,MD2_ATP); %return with ATP state of head 1
            
        else
            myosin = 3;
            life = -log(rand(1))/(kr2+kf4); %give the state a lifetime
            MD2(2) = -log(rand(1))/(kATP); %From state 1 so refresh ATP lifetime
            MD2(3) = MD2(3) + life; %add state3 lifetime to previous lifetime history
            MD1(4) = MD1(4) + life; %keep track of total time before ATP is released
            MD2(4) = MD2(4) + life; %keep track of total time before ATP is released
            myolife = myolife + life;
            [MD1,MD2,MD1_ATP,MD2_ATP] = ATP_rel(kATP,MD1,MD2,percent_ATP,MD1_ATP,MD2_ATP); %return with ATP state of head 2

        end
    end

    if myosin == 2
        myostate(2) = myostate(2) + 1;
        if rand(1)< (kf3/(kr1+kf3))
            myosin = 4;
            life = -log(rand(1))/(kr4+kr3); %give the state a lifetime
            MD2(2) = -log(rand(1))/(kATP); %From state 2 so refresh ATP lifetime
            MD2(3) = life; %new lifetime for MD2 because we've come from state 2
            MD1(3) = MD1(3) + life; %add state4 lifetime to previous lifetime history
            MD1(4) = MD1(4) + life; %keep track of total time before ATP is released
            MD2(4) = MD2(4) + life; %keep track of total time before ATP is released
            myolife = myolife + life; 
            [MD1,MD2,MD1_ATP,MD2_ATP] = ATP_rel(kATP,MD1,MD2,percent_ATP,MD1_ATP,MD2_ATP); %return with ATP state of heads 1&2
        else
            myosin = 1; 
        end
    end

   if myosin == 3
      myostate(3) = myostate(3) + 1;
      
        if rand(1)< (kf4/(kr2+kf4))
            myosin = 4;
            life = -log(rand(1))/(kr4+kr3); %give the state a lifetime 
            MD1(2) = -log(rand(1))/(kATP); %From state 3 so refresh ATP lifetime
            MD1(3) = life; %new lifetime for MD1 because we've come from state 3
            MD2(3) = MD2(3) + life; %add state4 lifetime to previous lifetime history
            MD1(4) = MD1(4) + life; %keep track of total time before ATP is released
            MD2(4) = MD2(4) + life; %keep track of total time before ATP is released
            myolife = myolife + life;
            [MD1,MD2,MD1_ATP,MD2_ATP] = ATP_rel(kATP,MD1,MD2,percent_ATP,MD1_ATP,MD2_ATP); %return with ATP state of heads 1&2
        else
            myosin = 1; 
        end
   end   

  if myosin == 4
     myostate(4) = myostate(4) + 1;
        if rand(1)< (kr3/(kr3+kr4))
            myosin = 2;
            life = -log(rand(1))/(kr1+kf3); %give the state a lifetime
            MD2(3) = 0; %reset md2 lifetime because it is down
            MD1(3) = MD1(3) + life; %add state2 lifetime to previous lifetime history
            MD1(4) = MD1(4) + life; %keep track of total time before ATP is released
            MD2(4) = MD2(4) + life; %keep track of total time before ATP is released
            myolife = myolife + life;
            [MD1,MD2,MD1_ATP,MD2_ATP] = ATP_rel(kATP,MD1,MD2,percent_ATP,MD1_ATP,MD2_ATP); %return with ATP state of head 1
        else
            myosin = 3; 
            life = -log(rand(1))/(kr2+kf4); %give the state a lifetime
            MD1(3) = 0; %reset md1 lifetime because it is down
            MD2(3) = MD2(3) + life; %add state3 lifetime to previous lifetime history
            MD1(4) = MD1(4) + life; %keep track of total time before ATP is released
            MD2(4) = MD2(4) + life; %keep track of total time before ATP is released
            myolife = myolife + life;
            [MD1,MD2,MD1_ATP,MD2_ATP] = ATP_rel(kATP,MD1,MD2,percent_ATP,MD1_ATP,MD2_ATP); %return with ATP state of head 2
        end
    end


end     

    end

