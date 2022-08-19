function [X, num_features, num_samples_total, num_clusters] = call_dataset()

load ovariancancer.mat
X = obs.';
[num_features,num_samples_total] = size(X);
num_clusters = 2;

end
