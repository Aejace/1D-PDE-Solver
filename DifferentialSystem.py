import math
import numpy
from scipy.integrate import solve_ivp
from BoundaryConditions import LeftFixedValueBoundaryCondition, LeftHeatFluxBoundaryCondition, \
    LeftCoolingBoundaryCondition, RightFixedValueBoundaryCondition, RightHeatFluxBoundaryCondition, \
    RightCoolingBoundaryCondition


class DifferentialSystem:
    # Define differential system parameters
    deltaX = 0.1
    deltaT = 0.1
    xSamplePoints = []
    timeSamplePoints = []
    initialState = []
    leftBoundaryCondition = None
    rightBoundaryCondition = None
    computationalSolution = None

    # TODO: Define PDE parameters
    alpha = 1.0
    lateralCoefficientOfCooling = 0.0
    lateralAmbientY = 5.0
    leftStartingY = 5.0
    rightStartingY = 30.0


    def __init__(self):
        self.SpecifySpaceParameters()
        self.SpecifyTimeParameters()
        self.SpecifyBoundaryConditions()

        # Future work: Could write user-gets for other parameters
        # self.SpecifyInitialDistributionFunction()
        # Specify other parameters ---

        self.SetInitialState()


    def SpecifySpaceParameters(self):
        # Future work: Could be parameterized, and implement scaling of the solution
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


    def SpecifyBoundaryConditions(self):
        # Print options
        print("Please select a left boundary condition type: ")
        print("1: Fixed value")
        print("2: Prescribed flux")
        print("3: Prescribed cooling")
        selection = input('Selection: ')

        # Validate
        while not (selection == "1" or selection == "2" or selection == "3"):
            selection = input('Error: Please input a number between 1 and 3: ')

        # Instantiate selected boundary condition
        leftBoundaryCondition = None
        match selection:
            case "1":
                self.leftBoundaryCondition = LeftFixedValueBoundaryCondition(self)

            case "2":
                self.leftBoundaryCondition = LeftHeatFluxBoundaryCondition(self)

            case "3":
                self.leftBoundaryCondition = LeftCoolingBoundaryCondition(self)

        # Print options
        print("")
        print("Please select a right boundary condition type: ")
        print("1: Fixed value")
        print("2: Prescribed flux")
        print("3: Prescribed cooling")
        selection = input('Selection: ')

        # Validate
        while not (selection == "1" or selection == "2" or selection == "3"):
            selection = input('Error: Please input a number between 1 and 3: ')

        # Instantiate selected boundary condition
        match selection:
            case "1":
                self.rightBoundaryCondition = RightFixedValueBoundaryCondition(self)

            case "2":
                self.rightBoundaryCondition = RightHeatFluxBoundaryCondition(self)

            case "3":
                self.rightBoundaryCondition = RightCoolingBoundaryCondition(self)


    def SetInitialState(self):
        # Instantiate empty array to hold initial distribution of Y values at start time
        n = len(self.xSamplePoints)
        self.initialState = numpy.zeros(n, dtype=float)

        # non-dimensional initial conditions
        # Future work: Could be replaced with a function parameter rather than using an explicit linear function here.
        for i in range(0, n):
            x = float(i) * self.deltaX
            initialValueAtPositionX = self.leftStartingY + ((self.rightStartingY - self.leftStartingY) * x)  # Linear starting Y values varying from left to right.
            self.initialState[i] = initialValueAtPositionX


    # Generates the system of ODEs that represent the 1D partial differential equation being solved
    def GenerateOrdinaryDifferentialEquationSystem(self, t, state):
        n = len(self.xSamplePoints)
        ODEs = numpy.zeros(n, dtype=float)
        simplifiedConstants = self.alphaOverDeltaXSquared()

        # TODO: Enter differential equations here to modify program
        for i in range(1, n - 1):
            ODEs[i] = simplifiedConstants * (state[i - 1] - 2.0 * state[i] + state[i + 1]) + self.lateralCoefficientOfCooling * (self.lateralAmbientY - state[i])

        # Set left boundary conditions
        ODEs[0] = self.leftBoundaryCondition.ODE(state)
        # Set right boundary conditions
        ODEs[n - 1] = self.rightBoundaryCondition.ODE(state)

        print('-> ', end='')

        return ODEs


    # Generates Jacobian matrix
    # Takes 3 arguments because the function is called in solve_ivp with 3 arguments. (Library requirements)
    def GenerateJacobian(self, t, q):  # implement the differential-equation jacobian in this function
        # Initialize empty n x n matrix
        n = len(self.xSamplePoints)
        jacobianMatrix = numpy.zeros((n, n), dtype=float)

        # Helper function to simplify algebra
        simplifiedConstants = self.alphaOverDeltaXSquared()

        # TODO: Enter partial derivatives of differential equations here to modify program
        for i in range(1, n - 1):
            jacobianMatrix[i][i - 1] = simplifiedConstants
            jacobianMatrix[i][i] = simplifiedConstants * (- 2.0) - self.lateralCoefficientOfCooling
            jacobianMatrix[i][i + 1] = simplifiedConstants

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


    # Define Helper Functions
    def alphaOverDeltaXSquared(self):
        return self.alpha / (self.deltaX * self.deltaX)

