% Specify the folder containing the simulations.
mainFolder = 

runFunctionInSubfolders(mainFolder);

function [] = get_snap(q,time)
%plot_bricks Plot of Brick Configuration
%   patch(X,Y,C)

ntime = size(q,2);          % number of time steps
nblocks = size(q,1)/3;      % number of blocks

w = 0.2*ones(nblocks,1);    % block width
% w(1) = 1;
h = 0.4*ones(nblocks,1);    % block height
% h(1) = .2;

X1 = zeros(ntime,nblocks);
X2 = zeros(ntime,nblocks);
Theta = zeros(ntime,nblocks);
for j = 1:ntime
    for i = 1:nblocks
        X1(j,i) = q(3*(i-1)+1,j);
        X2(j,i) = q(3*(i-1)+2,j);
        Theta(j,i) = q(3*i,j);
    end
end

box on
axis([-4*max(w) 4*max(w) -10*max(h)/2 (nblocks+1)*max(h)])
axis equal
xticks([])
yticks([])
hold on
for j = time
    axis([-4*max(w) 4*max(w) -10*max(h)/2 (nblocks+1)*max(h)])
    hold on
    for i = 1:nblocks
        x1 = X1(j,i)-w(i)/2*cos(Theta(j,i))+h(i)/2*sin(Theta(j,i));
        y1 = X2(j,i)-w(i)/2*sin(Theta(j,i))-h(i)/2*cos(Theta(j,i));
        x2 = X1(j,i)+w(i)/2*cos(Theta(j,i))+h(i)/2*sin(Theta(j,i));
        y2 = X2(j,i)+w(i)/2*sin(Theta(j,i))-h(i)/2*cos(Theta(j,i));
        x3 = X1(j,i)+w(i)/2*cos(Theta(j,i))-h(i)/2*sin(Theta(j,i));
        y3 = X2(j,i)+w(i)/2*sin(Theta(j,i))+h(i)/2*cos(Theta(j,i));
        x4 = X1(j,i)-w(i)/2*cos(Theta(j,i))-h(i)/2*sin(Theta(j,i));
        y4 = X2(j,i)-w(i)/2*sin(Theta(j,i))+h(i)/2*cos(Theta(j,i));
        block(i) = patch([x1 x2 x3 x4], [y1 y2 y3 y4],[0.9290 0.6940 0],...
            'EdgeColor',[0 0 0.5],'LineWidth',0.1);
    end
end

end

function [] = get_animation(q,movie_name)
%plot_bricks Plot of Brick Configuration
%   patch(X,Y,C)

ntime = size(q,2);          % number of time steps
nblocks = size(q,1)/3;      % number of blocks

w = 0.2*ones(nblocks,1);    % block width
% w(1) = 1;
h = 0.4*ones(nblocks,1);    % block height
% h(1) = .2;

Dtime = 0.0626418390534633;

X1 = zeros(ntime,nblocks);
X2 = zeros(ntime,nblocks);
Theta = zeros(ntime,nblocks);
for j = 1:ntime
    for i = 1:nblocks
        X1(j,i) = q(3*(i-1)+1,j);
        X2(j,i) = q(3*(i-1)+2,j);
        Theta(j,i) = q(3*i,j);
    end
end

v = VideoWriter(movie_name);
open(v);
figure(1)
clf
box on
% Y_max = max(X2(:,end))+1;
axis([-4*max(w) 4*max(w) -10*max(h)/2 (nblocks+1)*max(h)])
axis equal
hold on
set(gca,'XTick',[], 'YTick', [])
for j = 1:50:ntime
    text(3,5,num2str((j-1)*Dtime))
    axis([-4*max(w) 4*max(w) -10*max(h)/2 (nblocks+1)*max(h)])
    hold on
    for i = 1:nblocks
        x1 = X1(j,i)-w(i)/2*cos(Theta(j,i))+h(i)/2*sin(Theta(j,i));
        y1 = X2(j,i)-w(i)/2*sin(Theta(j,i))-h(i)/2*cos(Theta(j,i));
        x2 = X1(j,i)+w(i)/2*cos(Theta(j,i))+h(i)/2*sin(Theta(j,i));
        y2 = X2(j,i)+w(i)/2*sin(Theta(j,i))-h(i)/2*cos(Theta(j,i));
        x3 = X1(j,i)+w(i)/2*cos(Theta(j,i))-h(i)/2*sin(Theta(j,i));
        y3 = X2(j,i)+w(i)/2*sin(Theta(j,i))+h(i)/2*cos(Theta(j,i));
        x4 = X1(j,i)-w(i)/2*cos(Theta(j,i))-h(i)/2*sin(Theta(j,i));
        y4 = X2(j,i)-w(i)/2*sin(Theta(j,i))+h(i)/2*cos(Theta(j,i));
        block(i) = patch([x1 x2 x3 x4], [y1 y2 y3 y4],[0.9290 0.6940 0],...
            'EdgeColor',[0 0 0.5],'LineWidth',0.1);
    end
    frame = getframe(gcf);
    writeVideo(v,frame);
    pause(0.1)
    delete(block)
end

close(v);

end

function [] = get_bifurcation_analysis(fig_title)

    global h
    
    %% loading position vector of blocks
    load('q.mat')

    %% reshaping array
    n_bifurcations = size(q,1);
    if n_bifurcations>100
        n_bifurcations = 100;
    end
    n_dof = size(q,2);
    n_time = size(q,3);

    l = zeros(n_dof,n_time,n_bifurcations);

    for i = 1:size(q,1)
        for j = 1:size(q,2)
            for k = 1:size(q,3)
                l(j,k,i) = q(i,j,k);
            end
        end
        % get animation
        get_animation(l(:,:,i),['movie',int2str(i),'.avi'])
    end

    %% comparing motions of different blocks due to bifurcations
    figure()
    hold on
    box on
    hold on
    for i = 1:n_bifurcations
        title('$y$ coordinate of the CG of top block','Interpreter','latex')
        for j = n_dof-1
            if abs(l(n_dof-2,end,i)-l(n_dof-5,end,i))>0.2
                plot(l(j,:,i),'color','black')
            else
                plot(l(j,:,i),'color','black')
            end
        end
    end
    axis([1,size(l,2),l(n_dof-1,1,1)-0.2,l(n_dof-1,1,1)+0.2])
    saveas(gcf,[fig_title,'_y','.png'])
    
    %% comparing motions of different blocks due to bifurcations
    figure()
    hold on
    box on
    hold on
    for i = 1:n_bifurcations
        title('$\theta$ coordinate of the CG of top block','Interpreter','latex')
        for j = n_dof
            if abs(l(n_dof-2,end,i)-l(n_dof-5,end,i))>0.2
                plot(l(j,:,i),'color','black')
            else
                plot(l(j,:,i),'color','black')
            end
        end
    end
    saveas(gcf,[fig_title,'_theta','.png'])
    
    %% 
    figure(50)
    hold on
    box on
    hold on
    for i = 1:n_bifurcations
        title('$\theta$ coordinate of the CG of top block','Interpreter','latex')
        for j = n_dof
            plot(l(j,:,i),'linewidth',2)
        end
    end
    saveas(gcf,[fig_title,'_theta','.png'])

    %%
    figure()
    hold on
    box on
    title('$x$ coordinate of the CG of top block','Interpreter','latex')
    xlabel('iteration')
    for i = 1:n_bifurcations
        for j = n_dof-2
            if abs(l(n_dof-2,end,i)-l(n_dof-5,end,i))>0.2
                plot(l(j,:,i),'color','black')
            else
                plot(l(j,:,i),'color','black')
            end
        end
    end
    saveas(gcf,[fig_title,'_x','.png'])

    %% plot snapshot
    figure()
    hold on
    box on
    for i = 1:n_bifurcations
        subplot(ceil(sqrt(n_bifurcations)),ceil(sqrt(n_bifurcations)),i)
        box on
        hold on
        title(num2str(i))
        get_snap(l(:,:,i),size(l,2))
    end
    saveas(gcf,[fig_title,'_snap','.png'])
end

% Main script
function runFunctionInSubfolders(mainFolder)
    % List all folders and subfolders
    folders = dir(mainFolder);
    folders = folders([folders.isdir]);
    folders = folders(~ismember({folders.name}, {'.', '..'}));
    
    % Iterate over each subfolder
    theta_save = zeros(numel(folders),3);
    for i = 1:numel(folders)
        subfolder = fullfile(mainFolder, folders(i).name);
        cd(subfolder); % Change current directory to the subfolder
            fig_title = extractBetween(folders(i).name,"outputs","_mu");
            get_bifurcation_analysis(folders(i).name); % Run the function

    end
end

