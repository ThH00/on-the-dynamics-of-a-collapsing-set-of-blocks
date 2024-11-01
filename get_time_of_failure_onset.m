% specify folder constaining simulations
mainFolder = 
% Example usage
% mainFolder = 'Z:\chone\Desktop\theresa\stacked-blocks-2d\outputs\2024-06-04_13-43-21';

[fail_array] = runFunctionInSubfolders(mainFolder);
fail_array = sortrows(fail_array);
figure()
hold on
box on
plot(fail_array(:,1),1./fail_array(:,2),'linewidth',2);
xlabel('$\omega$','Interpreter','latex')
ylabel('$1/{text{number of iterations to fail}}$','Interpreter','latex')

figure()
hold on
dtime = 2*pi/100./fail_array(:,1);
fail_array(:,3) = fail_array(:,2).*dtime;
box on
plot(fail_array(:,1),1./fail_array(:,3),'linewidth',2);
xlabel('$\omega$','Interpreter','latex')
ylabel('$\frac{1}{t_F}$')


load('corners.mat')

%% functions

function [min_t] = locate_fail()
    load('corners.mat')

    n_bifurcations = size(corners,1);
    n_dof = size(corners,2);
    n_time = size(corners,3);

    t = zeros(n_bifurcations,1);
    for i = 1:n_bifurcations
        sample = corners(i,:,:);
        sample = reshape(sample,n_dof,n_time);
        % find first time same is 4
        [x,y] = find(sample==4,1);
        if size(x,1) == 0
            t(i) = n_time; %Inf;
        else
            t(i) = y;
        end
    end
    min_t = min(t);
end

function [fail_array] = runFunctionInSubfolders(mainFolder)
    % List all folders and subfolders
    folders = dir(mainFolder);
    folders = folders([folders.isdir]);
    folders = folders(~ismember({folders.name}, {'.', '..'}));
    
    % Iterate over each subfolder
    fail_array = zeros(numel(folders),2);
    for i = 1:numel(folders)
        subfolder = fullfile(mainFolder, folders(i).name);
        cd(subfolder); % Change current directory to the subfolder
        frequency = extractBetween(folders(i).name,"omega","_mu");
        frequency = cell2mat(frequency);
        frequency = str2num(frequency);
        t = locate_fail(); % Run the function
        fail_array(i,1) = frequency;
        fail_array(i,2) = t;
    end
end

