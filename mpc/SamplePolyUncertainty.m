function [output] = SamplePolyUncertainty(Vertices)
%SAMPLEPOLYUNCERTAINTY Summary of this function goes here
%   Vertices: 1 x N cell of vertices of a set
N = length(Vertices);
if N == 1
    output = Vertices{1};
    return
end

nx = size(Vertices{1}, 2);

lambda = sort(rand(1, N-1));
lambda = [0.0 lambda 1.0];
weights = diff(lambda);
output = cell2mat(Vertices)*kron(weights, eye(nx))';

end

