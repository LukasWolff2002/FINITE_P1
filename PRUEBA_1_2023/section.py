import numpy as np

class Section:
    def __init__(self, thickness, E, nu, type='planeStress'):
        self.thickness = thickness
        self.E = E
        self.nu = nu
        self.type = type
        self.D = self._compute_D()

    def _compute_D(self):
        E, nu = self.E, self.nu
        
        if isinstance(E, np.ndarray):
            return E
        
        if self.type == 'planeStress':
            return (E / (1 - nu**2)) * np.array([
                [1, nu, 0],
                [nu, 1, 0],
                [0, 0, (1 - nu) / 2]
            ])
        elif self.type == 'planeStrain':
            coef = E / ((1 + nu)*(1 - 2*nu))
            return coef * np.array([
                [1 - nu, nu, 0],
                [nu, 1 - nu, 0],
                [0, 0, (1 - 2*nu) / 2]
            ])
        else:
            raise ValueError(f"Invalid type: {self.type}")