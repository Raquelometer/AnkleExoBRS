function [c_j] = computeLevelCandidate(A_con, b_con, Qval)
%COMPUTELEVELCANDIDATE Given an elliposed x'Px and a constraint polytope, 
%  this function computes the level set  of the maximum circumscribed
%  ellipsoid

c_j = inf;
for j = 1:length(b_con)
    c_j = min(c_j,(b_con(j)^2 )/ (A_con(j,:)*Qval*A_con(j,:)'));
end
end