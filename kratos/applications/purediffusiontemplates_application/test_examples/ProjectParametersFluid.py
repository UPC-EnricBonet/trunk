domain_size = 2


# Fluid solver configuration
class SolverSettings:
    solver_type =  "fluid_solver"
    domain_size = 2
    #TurbulenceModel = "None"

    """# Velocity solver
    class velocity_linear_solver_config:
        solver_type = "BiConjugate gradient stabilized"
        tolerance = 1E-5
        max_iteration = 100
        preconditioner = "ILU0"
        scaling = False"""
    # Pressure solver
    """class pressure_linear_solver_config:
        solver_type = "AMGCL"
        tolerance = 1E-5
        max_iteration = 100
        preconditioner = "ILU0"
        scaling = False
        krylov_type = "CG"
        smoother_type = "DAMPED_JACOBI"
    
    #convergence criteria settings
    vel_toll = 1E-3
    press_toll = 1E-3
    divergence_cleareance_step = 10
    
    #other solver settings
    oss_switch = 1
    compute_reactions = True
    time_order = 2
    predictor_corrector = False
    dynamic_tau = 0.01
    max_vel_its = 4
    max_press_its = 3
    laplacian_form = 1"""

#general problem settings
AutomaticDeltaTime = "Fixed"
Dt = 0.05
Start_time = 0.0
max_time = 20
nsteps = 2000


groups_dictionary = {
        "area" : 1,
                   }
#output settings
output_time = 0.1
output_step = 10
VolumeOutput = True
nodal_results=["PRESSURE"]
gauss_points_results=[]
GiDPostMode = "Ascii"
GiDWriteMeshFlag = True
GiDWriteConditionsFlag = True
GiDWriteParticlesFlag = False
GiDMultiFileFlag = "Single"

problem_name="FluidTransport"
"""problem_path="/home/pbecker/ric/bilal.gid"
kratos_path="D:\Kratos"""
