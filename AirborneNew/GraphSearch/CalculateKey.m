function key = CalculateKey(node,g,rhs,goalNode,km)
% function key = CalculateKey(node,g,rhs,goalNode)    
key = [min([g(node),rhs(node)])+ h(node,goalNode)+ km , ...
       min([g(node),rhs(node)])];
end 