class LeftFixedValueBoundaryCondition:
    # Local Variables
    name = "Fixed value"
    PDE = None

    def __init__(self, PDE):
        self.PDE = PDE


    def Initialize(self):
        return


    def ODE(self, state):
        return 0


    def PartialDerivative(self, jacobianMatrix):
        n = len(jacobianMatrix[0])
        for i in range(0, n):
            # left boundary - fixed temperature
            jacobianMatrix[0][i] = 0.0
        return


class RightFixedValueBoundaryCondition:
    # Local Variables
    name = "Fixed value"
    PDE = None

    def __init__(self, PDE):
        self.PDE = PDE


    def Initialize(self):
        return


    def ODE(self, state):
        return 0


    def PartialDerivative(self, jacobianMatrix):
        n = len(jacobianMatrix[0])
        for i in range(0, n):
            # left boundary - fixed temperature
            jacobianMatrix[n - 1][i] = 0.0