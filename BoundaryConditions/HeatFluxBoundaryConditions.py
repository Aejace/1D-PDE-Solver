class LeftHeatFluxBoundaryCondition:
    # Local Variables
    name = "Heat flux"
    PDE = None
    leftFlux = None
    simplifiedConstants = None


    def __init__(self, PDE):
        self.PDE = PDE
        self.simplifiedConstants = PDE.alphaOverDeltaXSquared()


    def Initialize(self):
        self.leftFlux = float(input('Enter left heat flux value: ')) # Future work: Add error handling


    def ODE(self, state):
        return 2.0 * self.simplifiedConstants * (state[1] - state[0]) - (2.0 / self.PDE.deltaX) * self.leftFlux


    def PartialDerivative(self, jacobianMatrix):
        jacobianMatrix[0][0] = -2.0 * self.simplifiedConstants
        jacobianMatrix[0][1] = 2.0 * self.simplifiedConstants


class RightHeatFluxBoundaryCondition:
    # Local Variables
    name = "Heat flux"
    PDE = None
    rightFlux = None
    simplifiedConstants = None


    def __init__(self, PDE):
        self.PDE = PDE
        self.simplifiedConstants = PDE.alphaOverDeltaXSquared()


    def Initialize(self):
        self.rightFlux = float(input('Enter right heat flux value: ')) # Future work: Add error handling


    def ODE(self, state):
        n = len(state)
        return - 2.0 * self.simplifiedConstants * (state[n - 1] - state[n - 2]) - (2.0 / self.PDE.deltaX) * self.rightFlux


    def PartialDerivative(self, jacobianMatrix):
        n = len(jacobianMatrix[0])
        jacobianMatrix[n - 1][n - 2] = 2.0 * self.simplifiedConstants
        jacobianMatrix[n - 1][n - 1] = -2.0 * self.simplifiedConstants
