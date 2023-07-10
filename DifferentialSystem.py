import math
import numpy
from scipy.integrate import solve_ivp
import PartialDifferentialEquations
from InitialConditions.LinearInitialCondition import LinearInitialCondition
from PartialDifferentialEquations.BatemanBurgersPDE import BatemanBurgersPDE
from PartialDifferentialEquations.HeatPDE import HeatPDE


class DifferentialSystem:
    # Define differential system parameters
    PDE = None
    leftBoundaryCondition = None
    rightBoundaryCondition = None
    initialConditionFunction = None
    computationalSolution = None
    deltaX = 0.1
    deltaT = 0.1
    xSamplePoints = []
    timeSamplePoints = []
    initialState = []


    def __init__(self):
        self.SpecifySpaceParameters()
        self.SpecifyTimeParameters()
        self.SpecifyPDE()
        self.SpecifyBoundaryConditions()
        self.SpecifyInitialConditionFunction()
        self.SetInitialState()


    def SpecifySpaceParameters(self):
        # Future work: Could be parameterized, and implement scaling of the solution # Note: linear initial distribution relies on start and end being 0 and 1
        startingPosition = 0
        endingPosition = 1

        numberOfPoints = int(input('Enter quantity of evenly spaced sampling points: '))

        # Validate input
        while numberOfPoints < 2:
            numberOfPoints = int(input('Invalid response: Enter quantity of evenly spaced sampling points (Must be greater than 2): '))

        # Calculate deltaX (x-axis step size)
        self.deltaX = (endingPosition - startingPosition) / (float(numberOfPoints) - 1.0)

        # Instantiate array of sampling points along x-axis
        self.xSamplePoints = numpy.zeros(numberOfPoints)
        for i in range(0, numberOfPoints):
            self.xSamplePoints[i] = self.deltaX * i


    def SpecifyTimeParameters(self):
        startTime = float(input('Enter initial time: '))
        endTime = float(input('Enter final time: '))
        self.deltaT = float(input('Enter time step size DeltaT: ')) # Stored for use in solve_ivp library function call

        # Calculate number of time steps needed to cover specified time period.
        numberOfTimeSteps = ((endTime - startTime) / self.deltaT) + 1.0

        # If the specified time steps do not evenly divide the time period, the step size must be adjusted such that it does.
        numberOfTimeSteps = float(int(numberOfTimeSteps * 1000)/1000.0) # Truncate to 3 decimals using type conversion witchcraft
        roundedNumberOfTimeSteps = math.ceil(numberOfTimeSteps) # Round up to the nearest whole number
        if roundedNumberOfTimeSteps > numberOfTimeSteps:
            # Adjust time step
            self.deltaT = (endTime - startTime) / (roundedNumberOfTimeSteps - 1.0)
            print ("Warning: Time step size does not divide time interval evenly, and has been adjusted to: " + str(self.deltaT))

        # Instantiate array of sampling points along time-axis
        self.timeSamplePoints = numpy.zeros(roundedNumberOfTimeSteps)
        for i in range(0, roundedNumberOfTimeSteps):
            self.timeSamplePoints[i] = self.deltaT * i


    def SpecifyPDE(self):
        # Print options
        print("Please select a PDE: ")
        print("1: Heat")
        print("2: Bateman-Burgers")
        selection = input('Selection: ')

        # Validate
        while not (selection == "1" or selection == "2"):
            selection = input('Error: Please input a number listed above: ')

        # Instantiate selected boundary condition
        leftBoundaryCondition = None
        match selection:
            case "1":
                self.PDE = HeatPDE(self.deltaX)

            case "2":
                self.PDE = BatemanBurgersPDE(self.deltaX)


    def SpecifyBoundaryConditions(self):
        # Print left options
        print("Please select a left boundary condition type: ")
        leftBoundaryConditionOptions = self.PDE.GetLeftBoundaryConditions()
        for i in range(0, len(leftBoundaryConditionOptions)):
            print( str(i + 1) + ": " + leftBoundaryConditionOptions[i].name + ": ")
        selection = int(input('Selection: '))

        # Validate
        while selection < 1 or selection > len(leftBoundaryConditionOptions):
            selection = int(input('Error: Please input a number listed above: '))

        # Instantiate selected boundary condition
        self.leftBoundaryCondition = leftBoundaryConditionOptions[selection - 1]
        self.leftBoundaryCondition.Initialize()


        # Print right options
        print("Please select a right boundary condition type: ")
        rightBoundaryConditionOptions = self.PDE.GetRightBoundaryConditions(selection - 1)
        for i in range(0, len(rightBoundaryConditionOptions)):
            print(str(i + 1) + ": " + rightBoundaryConditionOptions[i].name + ": ")
        selection = int(input('Selection: '))

        # Validate
        while selection < 1 or selection > len(rightBoundaryConditionOptions):
            selection = int(input('Error: Please input a number listed above: '))

        # Instantiate selected boundary condition
        self.rightBoundaryCondition = rightBoundaryConditionOptions[selection - 1]
        self.rightBoundaryCondition.Initialize()


    def SpecifyInitialConditionFunction(self):
        # Print options
        print("Please select an initial condition function type: ")
        print("1: Linear")
        selection = input('Selection: ')

        # Validate
        while not (selection == "1"):
            selection = input('Error: Please input a number listed above: ')

        # Instantiate selected boundary condition
        match selection:
            case "1":
                self.initialConditionFunction = LinearInitialCondition(self)


    def SetInitialState(self):
        # Instantiate empty array to hold initial distribution of Y values at start time
        n = len(self.xSamplePoints)
        self.initialState = numpy.zeros(n, dtype=float)

        # non-dimensional initial conditions
        # Future work: Could be replaced with a function parameter rather than using an explicit linear function here.
        for i in range(0, n):
            x = float(i) * self.deltaX
            initialValueAtPositionX = self.initialConditionFunction.GetValue(x)
            self.initialState[i] = initialValueAtPositionX


    # Generates the system of ODEs that represent the 1D partial differential equation being solved
    def GenerateOrdinaryDifferentialEquationSystem(self, t, state):
        n = len(self.xSamplePoints)
        ODEs = numpy.zeros(n, dtype=float)

        # Fill ODEs
        for i in range(1, n - 1):
            ODEs[i] = self.PDE.ODE(state, i)

        # Set left boundary conditions
        ODEs[0] = self.leftBoundaryCondition.ODE(state)

        # Set right boundary conditions
        ODEs[n - 1] = self.rightBoundaryCondition.ODE(state)

        print('-> ', end='')

        return ODEs


    # Generates Jacobian matrix
    # Takes 3 arguments because the function is called in solve_ivp with 3 arguments. (Library requirements)
    def GenerateJacobian(self, t, state):  # implement the differential-equation jacobian in this function
        # Initialize empty n x n matrix
        n = len(self.xSamplePoints)
        jacobianMatrix = numpy.zeros((n, n), dtype=float)

        # Fill matrix
        for i in range(1, n - 1):
            self.PDE.PartialDerivative(jacobianMatrix, state, i)

        # Set left boundary conditions
        self.leftBoundaryCondition.PartialDerivative(jacobianMatrix)
        # Set right boundary conditions
        self.rightBoundaryCondition.PartialDerivative(jacobianMatrix)

        return jacobianMatrix


    def SolveSystem(self):
        # Format time range
        timeRange = numpy.zeros(2)
        timeRange[0] = self.timeSamplePoints[0]
        timeRange[1] = self.timeSamplePoints[len(self.timeSamplePoints) - 1]

        print('...working...')
        # Calculate solution with library call to solve_ivp
        self.computationalSolution = solve_ivp(self.GenerateOrdinaryDifferentialEquationSystem, timeRange, self.initialState, method ='Radau', t_eval = self.timeSamplePoints, dense_output = True, atol = 1.0e-9, rtol = 1.0e-9, jac = self.GenerateJacobian)
        return self.computationalSolution

