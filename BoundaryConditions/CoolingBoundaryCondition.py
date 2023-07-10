class LeftCoolingBoundaryCondition:
    # Local Variables
    name = "Cooling"
    PDE = None
    simplifiedConstants = None
    leftCoefficientOfCooling = None
    leftAmbientY = None


    def __init__(self, PDE):
        self.PDE = PDE
        self.simplifiedConstants = PDE.alphaOverDeltaXSquared()


    def Initialize(self):
        self.leftCoefficientOfCooling = float(input('Enter left coefficient of cooling value: ')) # Future work: Add error handling
        self.leftAmbientY = float(input('Enter temperature an infinite distance from left side of the sample space: '))  # Future work: Add error handling


    def ODE(self, state):
        return 2.0 * self.simplifiedConstants * (state[1] - state[0]) - (2.0 / self.PDE.deltaX) * self.leftCoefficientOfCooling * (state[0] - self.leftAmbientY)


    def PartialDerivative(self, jacobianMatrix):
        jacobianMatrix[0][0] = -2.0 * self.simplifiedConstants - (2.0 / self.PDE.deltaX) * self.leftCoefficientOfCooling
        jacobianMatrix[0][1] = 2.0 * self.simplifiedConstants


class RightCoolingBoundaryCondition:
    # Local Variables
    name = "Cooling"
    PDE = None
    simplifiedConstants = None
    rightCoefficientOfCooling = None
    rightAmbientY = None


    def __init__(self, PDE):
        self.PDE = PDE
        self.simplifiedConstants = PDE.alphaOverDeltaXSquared()


    def Initialize(self):
        self.rightCoefficientOfCooling = float(input('Enter left coefficient of cooling value: ')) # Future work: Add error handling
        self.rightAmbientY = float(input('Enter temperature an infinite distance from right side of the sample space: '))  # Future work: Add error handling


    def ODE(self, state):
        n = len(state)
        return - 2.0 * self.simplifiedConstants * (state[n - 1] - state[n - 2]) - (2.0 / self.PDE.deltaX) * self.rightCoefficientOfCooling * (state[n - 1] - self.rightAmbientY)


    def PartialDerivative(self, jacobianMatrix):
        n = len(jacobianMatrix[0])
        jacobianMatrix[n - 1][n - 2] = 2.0 * self.simplifiedConstants
        jacobianMatrix[n - 1][n - 1] = -2.0 * self.simplifiedConstants - (2.0 / self.PDE.deltaX) * self.rightCoefficientOfCooling