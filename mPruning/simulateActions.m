function an=simulateActions(pa)

an = sum(rand>[0 cumsum(pa)]);
