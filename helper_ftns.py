def in_birth_rate(r, K, Nt, N):
    return max(0, (r * (1 - (Nt/K)) * N))

def out_birth_rate(r, K, m, Nt, N):
    return max(0, (r * m * (1 - (Nt/K)) * N))