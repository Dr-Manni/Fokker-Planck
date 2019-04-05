clear
close all

%## Parameters### 
%Potential 1
k1=0.003;
a1=3.3;

%### Potential 2
k2=1;
a2=3.5;
%###
k12 =-50;
c=1;
%x1 : fast
%x2 : slow
V = @(x1, x2) k1 * (x1.^2 - a1^2).^4 + k2 * (x2.^2 - a2^2).^2 + k12 ./ sqrt((x1 - x2).^2 + c.^2);

xmin = -6;
xmax = 6;
ymin = xmin;
ymax = xmax;

x1 = xmin:0.1:xmax;
x2 = x1;
[X1, X2] = meshgrid(x1, x2);

figure
surf(X2, X1, min(V(X1, X2), 100))
shading interp
xlabel('x2')
ylabel('x1')

% Generate points
npoints = 3000;
VX2 = xmin + (xmax - xmin) * rand(npoints, 1);
VX1 = ymin + (ymax - ymin) * rand(npoints, 1);

figure
voronoi(VX2, VX1, 'k.')
hold on
plot(VX2, VX1, 'w.')
plot(VX2, VX1,  'ko','MarkerSize',0.6)
xlabel('x2')
ylabel('x1')
axis square
% Adjaceny matrix
X = [VX2, VX1];

[ v, c ] = voronoin (X);
A = sparse ( npoints, npoints );

for i = 1 : npoints
    for j = i + 1 : npoints
        s = size ( intersect ( c{i}, c{j} ) );
        
        if ( 1 < s(2) )
            A(i, j) = 1;
            A(j, i) = 1;
        end
        
    end
end

% Boltzmann weights
sigma = 15;
beta = 2 / sigma^2;
W = zeros(npoints,1);
for i = 1 : npoints
    W(i) = exp(-beta * V(VX1(i), VX2(i)));
end

Q = zeros(npoints);

for i = 1 : npoints
    for j = 1 : npoints
        Q(i,j) = sqrt(W(j) / W(i)) * A(i,j);
    end
end

Q(logical(eye(size(Q)))) = - sum(Q,2);

[evecsQ, evalsQ] = eigs(sparse(Q'), 5, 'lr');
evalsQ = diag(evalsQ);
evalsQ = -sort(-evalsQ);

% Plot eigenvectors
figure
for ev=1:4
    subplot(1,4,ev)
    cc = evecsQ(:,ev);
    cc = cc./sum(abs(cc));
    

    xlabel('x2')
    ylabel('x1')
    xlim([xmin xmax])
    ylim([xmin xmax])
    axis square
    
    for i = 1 : npoints
        patch(v(c{i},1), v(c{i},2), cc(i), 'edgecolor', 'none'); % use color i.
    end

    colormap(redblue)

    %set(gca, 'CLim', [-0.003, 0.003]);
    set(gca, 'CLim', [-max(max(abs(cc))), max(max(abs(cc)))]);
end
    
% FLUX
Qw = zeros(npoints);
for i=1:npoints
    for j=1:npoints
       
            Qw(i,j) = sqrt(exp(-beta * 1) / exp(-beta * 1)) * A(i,j);
       
    end
end

Qw(logical(eye(size(Qw)))) = - sum(Qw,2);


% ANALYTICAL ROW OF THE TRANSITION MATRIX
D = sqrt((VX2-VX2').^2 + (VX1-VX1').^2);

D = mean(mean(nonzeros(D.*full(A))));     
     
P = @(x, y, tau)  1/(2*pi*tau*sigma^2) * exp(- x.^2 / (2*sigma^2*tau) - y.^2 / (2*sigma^2*tau));

Pi = zeros(1, npoints);
tau = 20;
dt = 0.001;

for i=1:npoints
    Pi(i) = P(VX2(i), VX1(i), tau*dt)*D;
end

figure
xlabel('x2')
ylabel('x1')
xlim([xmin xmax])
ylim([xmin xmax])
axis square
    
for i = 1 : npoints
    patch(v(c{i},1), v(c{i},2), Pi(i), 'edgecolor', 'none'); % use color i.
end


[~,start] = min(sqrt((VX2 - 0).^2 + (VX1 - 0).^2));



%Flux diffusion
     

C = zeros(1, npoints);

sdt = sqrt(dt);

for rep = 1:1000000
    %%MSM

    eta = randn(tau-1, 2);

    x = zeros(tau, 2);


 



    for t = 2 : tau
        x(t, 1) = x(t-1, 1) + sigma*eta(t-1, 1)*sdt;
        x(t, 2) = x(t-1, 2) + sigma*eta(t-1, 2)*sdt;

    end
    %% Calculate row i of the matrix T


    [~, finish] = min(sqrt((VX2 - x(tau, 1)).^2 + (VX1 - x(tau, 2)).^2));
    
    C(1, finish) = C(1, finish) + 1;


end

Pi = C ./ sum(C);

figure
xlabel('x2')
ylabel('x1')
xlim([xmin xmax])
ylim([xmin xmax])
axis square
    
for i = 1 : npoints
    patch(v(c{i},1), v(c{i},2), Pi(i), 'edgecolor', 'none'); % use color i.
end

% Flux error

fluxes = linspace(1300,1400,100); %il flusso va da 101 a 399
err = zeros(length(fluxes),1);
u=0;
T = expm(Q*tau*dt);

for h=fluxes
    u=u+1;
    disp([num2str(h)])

    Th = T^h;

    err(u) = sqrt(sum((Pi  - Th(start,:)).^2)/npoints);
	%    err(u) = sqrt(sum((Pi(start) - Th(start,start)).^2)/npoints);

end

[m, idx] = min(err);
flux = fluxes(idx);

figure
plot(fluxes, err)

% ITS
taus = importdata('../fff0/data/tau.txt');
evals2 = exp(1.7*flux * taus * evalsQ(2) * dt);
evals3 = exp(1.7*flux * taus * evalsQ(3) * dt);
evals4 = exp(1.7*flux * taus * evalsQ(4) * dt);
evals5 = exp(1.7*flux * taus * evalsQ(5) * dt);


its2 = - taus ./ log(evals2)
its3 = - taus ./ log(evals3)
its4 = - taus ./ log(evals4)
its5 = - taus ./ log(evals5)