# FTS params
B = 0.1     # excluded volume used in FTS
E = 1.0    # Bjerrum length used in FTS
abar = 1.0  # smearing length used in FTS
Nref = 1.0  # reference N used in FTS
# MD simulation params
rcut = 6.5   # short-ranged distance cutoff
dt = 0.005  # timestep
kspace = 'pppm' # pppm or ewald solver
neighbin = 1.0 # neighbor list bin size
# for salts
salts = False  # Add salts? True or False
NCA = 0      # number of cations to add
NCC = 0      # number of anions to add
