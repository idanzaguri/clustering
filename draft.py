def log2(n):
    if n == 0:
        return 0
    return math.log2(n)
    
def entropy(P):
    assert sum(P) >= 1-1e-15 and sum(P) <= 1+1e-15, "Entropy: sum of all Pi("+str(sum(P))+") sould be 1!"
    return sum([ -p*log2(p) for p in P ])