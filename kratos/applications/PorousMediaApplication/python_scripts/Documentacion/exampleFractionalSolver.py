from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# -*- coding: utf-8 -*-
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.PorousMediaApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(HEAD_LEVEL)
    model_part.AddNodalSolutionStepVariable(CONCENTRATION)
    model_part.AddNodalSolutionStepVariable(CONCENTRATION_EXTERN)

    print("variables for the fractional step flowTrans solver added correctly")


def AddDofs(model_part):

    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(HEAD_LEVEL)
        node.AddDof(CONCENTRATION)

    print("dofs for the fractional step flowTrans solver added correctly")


class FrationalStepFlowTransSolver:
    #######################################################################
    def __init__(self,model_part,domain_size):
  
	self.model_part = model_part
	self.domain_size = domain_size
        
        # neighbour search
        #number_of_avg_elems = 10
        #number_of_avg_nodes = 10
        #self.neighbour_search = FindNodalNeighboursProcess(
            #model_part, number_of_avg_elems, number_of_avg_nodes)

        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

        #definition of the solvers
	self.basic_linear_solver = SkylineLUFactorizationSolver()

        #definition of the convergence criteria
	#tolerance for the solver
        self.conv_criteria = DisplacementCriteria(1e-5,1e-15)
        
    #######################################################################
    #######################################################################
    def Initialize(self):
        #(self.neighbour_search).Execute()
        #creating the solution strategy
        CalculateReactionFlag = False
        ReformDofSetAtEachStep = False	
        MoveMeshFlag = False

        import strategyBasicFlow
        self.solver = strategyBasicFlow.SolvingStrategyPython(self.model_part,self.time_scheme,self.basic_linear_solver,self.conv_criteria,CalculateReactionFlag,ReformDofSetAtEachStep,MoveMeshFlag)
	self.solver.SetEchoLevel(3)

        self.solver.Check()

        #self.solver.ApplyFractionalVelocityFixity()

        #######################################################################
        #######################################################################

    def Solve(self):
        #if(self.ReformDofAtEachIteration):
            #self.solver.ApplyFractionalVelocityFixity()
        
        t1 = timer.time()
####################################################################################        
        

        # flow equation iterations
        number_of_flow_iterations = 1
        for i in range(0, number_of_pressure_iterations):
            print("pressure iteration number: ", i + 1)
            self.CalculateFlowEquationIteration(i + 1)
        t2 = timer.time()


        number_of_transport_iterations = 2
        for i in range(1, number_of_transport_iterations):
            print("pressure iteration number: ", i + 1)
            self.CalculateTransportEquationIteration(i + 1)

         t3 = timer.time()

#######################################################################################
 def CalculateFlowEquationIteration(self, iteration_number):
        self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 1)
        self.model_part.ProcessInfo.SetValue(NL_ITERATION_NUMBER, iteration_number)
        #if (iteration_number == 1):
            #(self.VariableUtils).CopyVectorVar(MESH_VELOCITY, VELOCITY, self.model_part.Nodes)
        (self.solver).Solve()  # implicit resolution of the pressure system. All the other tasks are explicit

        # setting the Fract step number to 2 to calculate the head level gradient effect on the darcy flow
        #self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 2)
        #(self.ExplicitStrategy).InitializeSolutionStep(self.model_part.ProcessInfo);
        #(self.ExplicitStrategy).AssembleLoop(self.model_part.ProcessInfo);
        #(self.ExplicitStrategy).FinalizeSolutionStep(self.model_part.ProcessInfo);

 def CalculateTransportEquationIteration(self, iteration_number):
        self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 2)
        self.model_part.ProcessInfo.SetValue(NL_ITERATION_NUMBER, iteration_number)
        #if (iteration_number == 1):
            #(self.VariableUtils).CopyVectorVar(MESH_VELOCITY, VELOCITY, self.model_part.Nodes)
        (self.solver).Solve()  # implicit resolution of the pressure system. All the other tasks are explicit


    # Para calcular otras variable dependiente de flujo o transporte sería así:
    #def CalculateExplicitViscosityContribution(self):
        #self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 0)  # explicit contribution by viscosity is defined as fract step = 0
        #(self.ExplicitStrategy).InitializeSolutionStep(self.model_part.ProcessInfo);
        #(self.ExplicitStrategy).AssembleLoop(self.model_part.ProcessInfo);
        #(self.ExplicitStrategy).FinalizeSolutionStep(self.model_part.ProcessInfo);

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)
########################################################################################################3
