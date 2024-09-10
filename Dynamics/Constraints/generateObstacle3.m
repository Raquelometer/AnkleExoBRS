function [obsFuncMatrix] = generateObstacle3(constraints, vs, axis_dims)
% Function for generating an obstacle given an array of constraints
% Inputs: CONSTRAINTS an array of function handles that express the constraints, of the
% form c(x) < 0. VS a cell array of vectors of states at which the constraints should be
% evaluated.

% obsFuncMatrix should have same dimensions as value function.

num_constraints = length(constraints);

%obsFuncMatrix = zeros(axis_dims);
obsFuncMatrix = ones(axis_dims)*Inf;

for i = 1:num_constraints
    for j = 1:axis_dims(1)
        for k = 1:axis_dims(2)
            for l = 1:axis_dims(3)
                if i == 1
                    obsFuncMatrix(j,k,l) = constraints{i}(vs{1}(j),vs{2}(k),vs{3}(l));
                else
                    obsFuncMatrix(j,k,l) = min(obsFuncMatrix(j,k,l), constraints{i}(vs{1}(j),vs{2}(k),vs{3}(l)));
                end
            end
        end

    end

end

    
end