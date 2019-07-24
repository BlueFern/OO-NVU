
load two_sec_stim_data.mat

% ****_air are experiments when the animal is breathing air
% ***_oxy is when the animal is breathing oxygen 
% ***_Fract - is the changes turned into fractional - without Fract they are Micromolar changes

% 19 = is the number of sessions (this data is from 5 animals the first four have 4 sessions the last animal has 3 - this is all in order (Animal 1=1:4, Animal 2=5 to 8) - There is one month gap between each of the sessions.
% 4 = regions from where the responses are from (1 = whisker region, 2 = Artery, 3 = Vein, 4 = Parenchyma)
% 200 = data length (time)
% 30 = the number of trials in each of the 19 sessions

% Regions
regions = {'Whisker','Artery','Vein','Parenchyma'};

test = zeros(size(tim));
sum_all = zeros(size(tim));

%%%%%% What do you want to plot? %%%%%%%%
plot_variable = Hbt_2s_oxy_fract;

% Mouse 1
sum_mouse = zeros(size(tim));
for region = 1:4
    sum_region = zeros(size(tim));
    
    for session = 1:4
        for trial = 1:30
            for time = 1:200
                test(time) = plot_variable(session,region,time,trial);
                sum_region(time) = sum_region(time) + test(time); % Sum up over region
                sum_all(time) = sum_all(time) + test(time); % Sum up everything
                sum_mouse(time) = sum_mouse(time) + test(time); % Sum up over mouse
            end
          
          % Plot every trial on one figure, goes through all mice
            figure(1000)
            hold all
            plot(tim, test)
            title('All mice/trials/regions')
       end
    end
    
    mean_region = sum_region/(30*4);
    
    figure(10)
    subplot(2,2,region)
    title(regions(region))
    hold all
    plot(tim, mean_region)
end
mean_mouse = sum_mouse/(30*4*4);
figure(11)
subplot(2,3,1)
hold all
plot(tim, mean_mouse)
title('Mouse 1')

% Mouse 2
sum_mouse = zeros(size(tim));
for region = 1:4
    sum_region = zeros(size(tim));
    
    for session = 5:8
        for trial = 1:30
            for time = 1:200
                test(time) = plot_variable(session,region,time,trial);
                sum_region(time) = sum_region(time) + test(time); % Sum up over region
                sum_all(time) = sum_all(time) + test(time); % Sum up everything
                sum_mouse(time) = sum_mouse(time) + test(time); % Sum up over mouse
            end
          
          % Plot every trial on one figure
            figure(1000)
            hold all
            plot(tim, test)
       end
    end
    
    mean_region = sum_region/(30*4);
    
    figure(10)
    subplot(2,2,region)
    title(regions(region))
    hold all
    plot(tim, mean_region)
end
mean_mouse = sum_mouse/(30*4*4);
figure(11)
subplot(2,3,2)
hold all
plot(tim, mean_mouse)
title('Mouse 2')

% Mouse 3
sum_mouse = zeros(size(tim));
for region = 1:4
    sum_region = zeros(size(tim));
    
    for session = 9:12
        for trial = 1:30
            for time = 1:200
                test(time) = plot_variable(session,region,time,trial);
                sum_region(time) = sum_region(time) + test(time); % Sum up over region
                sum_all(time) = sum_all(time) + test(time); % Sum up everything
                sum_mouse(time) = sum_mouse(time) + test(time); % Sum up over mouse
            end
          
          % Plot every trial on one figure
            figure(1000)
            hold all
            plot(tim, test)
       end
    end
    
    mean_region = sum_region/(30*4);
    
    figure(10)
    subplot(2,2,region)
    title(regions(region))
    hold all
    plot(tim, mean_region)
end
mean_mouse = sum_mouse/(30*4*4);
figure(11)
subplot(2,3,3)
hold all
plot(tim, mean_mouse)
title('Mouse 3')

% Mouse 4
sum_mouse = zeros(size(tim));
for region = 1:4
    sum_region = zeros(size(tim));
    
    for session = 13:16
        for trial = 1:30
            for time = 1:200
                test(time) = plot_variable(session,region,time,trial);
                sum_region(time) = sum_region(time) + test(time); % Sum up over region
                sum_all(time) = sum_all(time) + test(time); % Sum up everything
                sum_mouse(time) = sum_mouse(time) + test(time); % Sum up over mouse
            end
          
          % Plot every trial on one figure
            figure(1000)
            hold all
            plot(tim, test)
       end
    end
    
    mean_region = sum_region/(30*4);
    
    figure(10)
    subplot(2,2,region)
    title(regions(region))
    hold all
    plot(tim, mean_region)
end
mean_mouse = sum_mouse/(30*4*4);
figure(11)
subplot(2,3,4)
hold all
plot(tim, mean_mouse)
title('Mouse 4')

% Mouse 5
sum_mouse = zeros(size(tim));
for region = 1:4
    sum_region = zeros(size(tim));
    
    for session = 17:19
        for trial = 1:30
            for time = 1:200
                test(time) = plot_variable(session,region,time,trial);
                sum_region(time) = sum_region(time) + test(time); % Sum up over region
                sum_all(time) = sum_all(time) + test(time); % Sum up everything
                sum_mouse(time) = sum_mouse(time) + test(time); % Sum up over mouse
            end
          
          % Plot every trial on one figure
            figure(1000)
            hold all
            plot(tim, test)
       end
    end
    
    mean_region = sum_region/(30*3);
    
    figure(10)
    subplot(2,2,region)
    title(regions(region))
    hold all
    plot(tim, mean_region)
end
mean_mouse = sum_mouse/(30*4*3);
figure(11)
subplot(2,3,5)
hold all
plot(tim, mean_mouse)
title('Mouse 5')

% Plot mean of all animals and trials and regions
mean_all = sum_all/(30*4*19);
figure(12)
plot(tim, mean_all)
title('Average all mice/trials/regions')

% Save variable (depends on which 4D Hb array used! Change the variable name)
HbT_oxy = mean_all;

% Clear excess variables
clear trial time sum_mouse sum_region session region mean_mouse mean_region test sum_all mean_all