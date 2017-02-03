function h = h(node1, node2)
% function h = h(node1, node2)
% Description: h is an underestimate (=<) of the cost function c(node1,
% node2). h stands for "heuristics"
% For this example we take a simple norm_1:
h = abs(node2-node1);
end 