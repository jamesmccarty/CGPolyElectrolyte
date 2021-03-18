# parameters
vex = 0.006804     # excluded volume (v/b^3)
lB = 0.33    # Bjerrum length (lB/b)
atilde = 0.408248  # smearing length (a/b)
Nref = 1.0  # reference N used in FTS
rho =  0.8889     # density (only used to conver to FTS (rho*b^3))
scale_Rg = False    # If True, all distances will be scaled by Rg = Sqrt(Nref)/Sqrt(6)
# MD simulation params
rcut = 6.5   # short-ranged distance cutoff
dt = 0.005  # timestep
kspace = 'pppm' # pppm or ewald solver
neighbin = 1.0 # neighbor list bin size
# for salts
salts = False  # Add salts? True or False
NCA = 0      # number of cations to add
NCC = 0      # number of anions to add
