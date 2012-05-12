hc = bf_scheme('UD', [1,2], [1,0;0,1]);
hs = bf_scheme('SIR_Kalman', [1,2], [1,0;0,1], 1000, 1);

[x,X] = bf_cov(hc)
[x,X] = bf_cov(hs)
figure; hold on; axis equal
S = bf_sample(hs); plot(S(:,1),S(:,2),'g.')

bf_observe_linear_uncorrelated(hc, [1,0;0,1], [2,2], [0.1,0.1]);
bf_observe_linear_uncorrelated(hs, [1,0;0,1], [2,2], [0.1,0.1]);
[x,X] = bf_cov(hc)
[x,X] = bf_cov(hs)
S = bf_sample(hs); plot(S(:,1),S(:,2),'r.')

% Sampling results for hs
bf_unique_samples(hs)
bf_stochastic_samples(hs)

bf_predict_linear(hc,[2,1;0,1],[1,0;0,1],[1,0.1]);
bf_predict_linear(hs,[2,1;0,1],[1,0;0,1],[1,0.1]);
%bf_predict_linrz(hc, 'f_double', [2,1;0,1],[1,0;0,1],[1,0.1]);
%bf_predict_linrz(hs, 'f_double', [2,1;0,1],[1,0;0,1],[1,0.1]);
[x,X] = bf_cov(hc)
[x,X] = bf_cov(hs)
S = bf_sample(hs); plot(S(:,1),S(:,2),'k.')

%bf_observe_likelihood(hs, 'Lgauss', [2,2]);
%S = bf_sample(hs); plot(S(:,1),S(:,2),'r.')
