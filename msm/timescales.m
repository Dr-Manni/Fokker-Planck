clear 
close all

model='fff0';
data='data'
taus=importdata([model, '/', data, '/', 'tau.txt']);
taus = taus';
EVALS = zeros(3, length(taus));

k = 0;
for tau = taus
	k = k+1
	C = importdata([model, '/', data, '/Cxy', num2str(tau), '.txt']);
	
	m = length(C);
	lm = sqrt(m);
	T=zeros(m);

	s = sum(C, 2);
	valid = s > 0;
	T( valid,: ) = bsxfun(@rdivide, C(valid,:), s(valid) );

	%colormap(cmap)
    


	[~, evals]=eigs(T', 4, 'lr');
	evals = diag(real(evals));
	evals = -sort(-evals);
	
	EVALS(:,k) = evals(2:4);
end
	
	
	
	figure
	subplot(1,2,1)
	plot(taus, EVALS, 'o-', 'linewidth', 2)
	xlabel('lag-time (timesteps)')
	ylabel('Eigenvalue')
	legend('2nd eval', '3rd eval', '4th eval')
	
	ITS = - repmat(taus, [3,1])./log(EVALS);
	
	subplot(1,2,2)
	plot(taus, ITS, 'linewidth', 2)
	hold on
	plot(taus, taus, 'k')

	xlabel('lag-time (timesteps)')
	ylabel('ITS')