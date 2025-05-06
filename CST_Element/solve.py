import numpy as np

class Solve:
    def __init__(self, nodes, elements):
        self.nodes = nodes

           
        self.elements = elements
        self.ndof = len(nodes) * 2   # 2 DOFs por nodo
        self.K_global = np.zeros((self.ndof+1, self.ndof+1))
        self.f_global = np.zeros((self.ndof+1, 1))
        self.u_global = np.zeros((self.ndof+1, 1))

    def assemble(self):
        for elem in self.elements:
            ke = elem.Ke
            idx = elem.calculate_indices()
            for i in range(6):
                for j in range(6):
                    self.K_global[idx[i], idx[j]] += ke[i, j]

    def apply_force(self, dof_index, value):
        self.f_global[dof_index] += value

    def apply_forces_vector(self, force_vector):
        self.f_global += force_vector.reshape(-1, 1)

    def apply_boundary_conditions(self):
        self.fixed_dofs = []
        self.free_dofs = []

        for node in self.nodes:
            for dof_val, dof_idx in zip(node.restrain, node.dofs):
                if dof_val == 1:
                    self.fixed_dofs.append(dof_idx)
                else:
                    self.free_dofs.append(dof_idx)

        self.fixed_dofs = np.array(self.fixed_dofs)
        self.free_dofs = np.array(self.free_dofs)

        #print("N° de DOFs libres:", len(self.free_dofs))
        #print("N° de DOFs fijos:", len(self.fixed_dofs))


        # Aplicar condiciones
        for dof in self.fixed_dofs:
            self.K_global[dof, :] = 0
            self.K_global[:, dof] = 0
            self.K_global[dof, dof] = 1
            self.f_global[dof] = 0

    def check_zero_rows(self):
        zero_rows = np.where(~self.K_global.any(axis=1))[0]
        if len(zero_rows) > 0:
            pass
            #print(f"❌ Fila(s) completamente nulas en K_global: {zero_rows}")
        else:
            pass
            #print("✅ No hay filas nulas en K_global")

        return zero_rows


    def solve(self):

        self.assemble()

        # Guardar copias originales para calcular reacciones después
        self.K_original = self.K_global.copy()
        self.f_original = self.f_global.copy()

        self.apply_boundary_conditions()

        # Detectar DOFs realmente usados
        used_dofs = sorted(set(dof for node in self.nodes for dof in node.dofs))

        # Submatrices compactas
        K_reduced = self.K_global[np.ix_(used_dofs, used_dofs)]
        f_reduced = self.f_global[used_dofs]

        # Resolver sistema reducido
        u_reduced = np.linalg.solve(K_reduced, f_reduced)

        # Armar solución global
        self.u_global = np.zeros_like(self.f_global)
        self.u_global[used_dofs] = u_reduced
        return self.u_global

    def get_displacement_at_node(self, node_id):
        node = self.nodes[node_id - 1]  # Asume que el ID parte en 1
        ux = self.u_global[node.dofs[0], 0]
        uy = self.u_global[node.dofs[1], 0]
        return ux, uy

    def print_summary(self):
       
        for i, node in enumerate(self.nodes):
            ux, uy = self.get_displacement_at_node(node.id)
           

    def compute_reactions(self):

        # Usar las matrices antes de aplicar las condiciones de borde
        R_total = self.K_original @ self.u_global - self.f_original

        self.reactions = np.zeros_like(self.f_global)
        self.reactions[self.fixed_dofs] = R_total[self.fixed_dofs]

        return self.reactions
    
    def print_applied_forces(self):

        # Usa f_original si existe, si no, usa f_global
        f = self.f_original if hasattr(self, 'f_original') else self.f_global

        for node in self.nodes:
            dof_x, dof_y = node.dofs
            fx = f[dof_x][0] if dof_x < len(f) else 0.0
            fy = f[dof_y][0] if dof_y < len(f) else 0.0

    

            






