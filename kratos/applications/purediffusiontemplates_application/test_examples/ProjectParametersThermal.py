domain_size = 2

Linear = True

Stationary=True

class SolverSettings1:
    solver_type = "thermal_solver"
    domain_size= 2
    time_order = 2
    predictor_corrector = False
    ReformDofAtEachIteration = False
    echo_level=0


    ###set variables to be used
    unknown_variable = "CONCENTRATION"
    density_variable= "DENSITY"
    """projection_variable= "TEMP_CONV_PROJ"
    volume_source_variable= "CONCENTRATION_FLUX"
    diffusion_variable= "CONDUCTIVITY"
    surface_source_variable= "FACE_CONCENTRATION_FLUX"
    mesh_velocity_variable= "MESH_VELOCITY"
    velocity_variable= "VELOCITY"""
    #specific_heat_variable= "SPECIFIC_HEAT"


    # Temperature solver
    class linear_solver_config:
        solver_type = "BiConjugate gradient stabilized"
        tolerance = 1E-6
        max_iteration = 5000
        preconditioner = "ILU0"
        scaling = False


"""class SolverSettings2:
    solver_type = "nonlinear_convection_diffusion_solver"
    domain_size= 2
    time_order = 2
    predictor_corrector = False
    ReformDofAtEachIteration = False
    echo_level=0
    
    #convergence criteria settings
    max_iter = 15
    toll = 1E-3


    ###set variables to be used
    unknown_variable = "TEMPERATURE"
    density_variable= "DENSITY"
    projection_variable= "TEMP_CONV_PROJ"
    volume_source_variable= "HEAT_FLUX"
    diffusion_variable= "CONDUCTIVITY"
    surface_source_variable= "FACE_HEAT_FLUX"
    mesh_velocity_variable= "MESH_VELOCITY"
    velocity_variable= "VELOCITY"
    specific_heat_variable= "SPECIFIC_HEAT"


    # Temperature solver
    class linear_solver_config:
        solver_type = "BiConjugate gradient stabilized"
        tolerance = 1E-6
        max_iteration = 5000
        preconditioner = "ILU0"
        scaling = False"""

nodal_results=["CONCENTRATION"]
gauss_points_results=[]
GiDPostMode = "Ascii"
GiDWriteMeshFlag = True
GiDWriteConditionsFlag = False
GiDMultiFileFlag = "Single"

AutomaticDeltaTime = "Fixed"
Dt = 0.05
Start_time = 0.0
max_time = 1.0
nsteps = 1000

#output settings
output_time = 0.1
output_step = 100
VolumeOutput = True

"""problem_name="bilal_convdiffConvectionDiffusion"
problem_path="/home/pbecker/ric/bilal_convdiff.gid"
kratos_path="D:\Kratos"""
