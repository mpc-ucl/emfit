function [a] = simulateEffort(pa); 

	a = sum(rand>[0 cumsum(pa')]);
