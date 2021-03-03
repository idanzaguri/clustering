def log2(n):
    if n == 0:
        return 0
    return math.log2(n)
    
def entropy(P):
    assert sum(P) >= 1-1e-15 and sum(P) <= 1+1e-15, "Entropy: sum of all Pi("+str(sum(P))+") sould be 1!"
    return sum([ -p*log2(p) for p in P ])
    
    
"""
print(type(idan_fg_clustering))
pk_file = open('c1.pk','wb')
pickle.dump(idan_fg_clustering, pk_file)
pk_file.close()
"""

pk_file = open('c1.pk','rb')
raw_data = pickle.load(pk_file)
pk_file.close()
print(type(raw_data))
