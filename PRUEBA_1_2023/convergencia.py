#===================================
#LOGICA DE LA IMPLEMENTACION
#===================================

# 1. En primer lugar es nesesario crear la malla sobre el archivo .geo
#    para eso, usar la funcion generate_mesh

# 2. Luego, es nesesario crear los grupos de nodos, segun se define en 
#    el archivo geo, de este modo, es posible definir las restricciones y cargas

# 3. Una vez definidos los nodos, es posible definir las secciones
#    donde el espesor de la seccion dependera del tipo de grupo

#===================================
#IMPORTANTE A TENER EN CUENTA
#===================================

# 1. La numeracion de nodos y elementos debe partir desde 1

#====================================
#FUNCIONES EXISTENTES
#====================================

# 1. Node, se encarga de crear los nodos
#    Node(id, x, y, restrain)
#    restrain por default es libre, pero se puede definir como
#    restrain = [i, j]

# 2. Section, se encarga de crear las secciones
#    Section(thickness, E, nu)
#    por default, section es planeStress, pero se puede definir como
#    section = Section(thickness, E, nu, type='planeStrain')

# 3. CST, se encarga de crear los elementos CST
#    CST(id, [node1, node2, node3], section)
#    Donde los nodos correspondientes de cada elemento se obtiene a
#    partir de la malla generada en el archivo .msh

# 4. Solve, se encarga de resolver el sistema, en este caso, en la primera instancia se hace
#    estructure = Solve(nodes, elements)
#    Donde nodes es la lista de nodos y elements es la lista de elementos

# 5. Graph es un archivo que contiene todas las funciones necesarias para graficar

#====================================
#IMPORTO LAS DEPENDENCIAS
#====================================

#Librerias base
import os
import numpy as np
import matplotlib.pyplot as plt
import gmsh
import meshio
import time

#Librerias creadas
from nodes import Node
from section import Section
from cst import CST
from solve import Solve
from graph import plot_all_elements, plot_applied_forces, plot_deformed_structure, plot_deformed_with_reactions, plot_von_mises_field, plot_all_scalar_fields_separately, plot_principal_fields,plot_von_mises_per_element, plot_elements_by_thickness, plot_results, compute_nodal_von_mises

#====================================
#FUNCIONES
#====================================

def generate_mesh(input_file, output_file, lc=5):
    gmsh.initialize()
    gmsh.open(input_file)

    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", lc)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)
    gmsh.write(output_file)
    gmsh.finalize()

def make_nodes_groups(output_file, title, restrxy=None, restrx=None, restry=None):
    mesh = meshio.read(output_file)
    tag_to_name = {v[0]: k for k, v in mesh.field_data.items()}
    grupos = {}

    # Elementos tipo "triangle"
    for cell_block, phys_tags in zip(mesh.cells, mesh.cell_data["gmsh:physical"]):
        if cell_block.type != "triangle":
            continue
        for tri, tag in zip(cell_block.data, phys_tags):
            nombre = tag_to_name.get(tag, f"{tag}")
            if nombre not in grupos:
                grupos[nombre] = []
            for node_id in tri:
                x, y = mesh.points[node_id][:2]
                grupos[nombre].append(Node(node_id + 1, x, y))

    # Elementos tipo "line" (por ejemplo, para restricciones o cargas)
    for cell_block, phys_tags in zip(mesh.cells, mesh.cell_data["gmsh:physical"]):
        if cell_block.type != "line":
            continue
        for line, tag in zip(cell_block.data, phys_tags):
            nombre = tag_to_name.get(tag, f"{tag}")
            if nombre not in grupos:
                grupos[nombre] = []
            for node_id in line:
                x, y = mesh.points[node_id][:2]
                restrain = [0, 0]
                if restrxy is not None and nombre in restrxy:
                    restrain = [1, 1]
                elif restrx is not None and nombre in restrx:
                    restrain = [1, 0]
                elif restry is not None and nombre in restry:
                    restrain = [0, 1]
                grupos[nombre].append(Node(node_id + 1, x, y, restrain=restrain))

    # Eliminar duplicados dentro de cada grupo
    for nombre in grupos:
        nodos_unicos = {}
        for nodo in grupos[nombre]:
            nodos_unicos[nodo.id] = nodo
        grupos[nombre] = list(nodos_unicos.values())

    # Visualizar nodos por grupo
    Node.plot_nodes_por_grupo(grupos, title, show_ids=False)

    return grupos, mesh

def make_sections(grupos, thickness, E, nu):
    sections = {}

    # Propiedades del material (ortotrópico PLA impreso)
    for group in thickness:
        sections[group] = Section(thickness[group], E, nu)

    # Diccionario global de nodos para búsqueda por ID
    nodes_dict = {}
    for group in grupos:
        for node in grupos[group]:
            nodes_dict[node.id] = node

    return sections, nodes_dict

def make_cst_elements (mesh, sections, nodes_dict):

    triangles = mesh.cells_dict['triangle']  # nodos por cada triángulo
    tags = mesh.cell_data_dict["gmsh:physical"]["triangle"]
    elements = []
    nodes = set()  # Usamos un set para evitar nodos repetidos

    for i in range(len(tags)):
        section = sections[str(tags[i])]

        node_ids = triangles[i]
        nodo_a = nodes_dict[node_ids[0]+1]
        nodo_b = nodes_dict[node_ids[1]+1]
        nodo_c = nodes_dict[node_ids[2]+1]

        for nodo in [nodo_a, nodo_b, nodo_c]:
            nodes.add(nodo)  # Se agregan al set automáticamente sin duplicados

        elem = CST(i+1, [nodo_a, nodo_b, nodo_c], section)
        elements.append(elem)

    nodes = list(nodes)

    return elements, nodes

def apply_distributed_force_y(grupo_nodos, fuerza_total_y, estructura):

    # Paso 1: ordena nodos si es necesario (aquí asumimos que ya están ordenados)
    nodos = grupo_nodos
    n = len(nodos)
    if n < 2:
        print("Se requieren al menos dos nodos para aplicar fuerza distribuida.")
        return

    # Paso 2: calcular longitud total de la línea
    longitudes = []
    total_length = 0
    for i in range(n - 1):
        dx = nodos[i+1].x - nodos[i].x
        dy = nodos[i+1].y - nodos[i].y
        L = np.sqrt(dx**2 + dy**2)
        longitudes.append(L)
        total_length += L

    # Paso 3: calcular carga distribuida por unidad de longitud
    q_total = fuerza_total_y  # Fuerza total a repartir
    q_lineal = q_total / total_length  # N/m

    # Paso 4: aplicar cargas parciales a cada nodo (2 nodos por segmento)
    nodal_forces = {node.id: np.array([0.0, 0.0]) for node in nodos}

    for i in range(n - 1):
        ni = nodos[i]
        nj = nodos[i + 1]
        xi, yi = ni.x, ni.y
        xj, yj = nj.x, nj.y

        dx = xj - xi
        dy = yj - yi
        L = longitudes[i]

        if abs(dx) < 1e-8:  # Tramo casi vertical
            # Aplica carga normal (hacia la izquierda en este caso)
            fx = 0.0
            fy = -q_lineal * L
        else:
            # Dirección normalizada del tramo
            vx = dx / L
            vy = dy / L

            # Vector normal unitario (perpendicular al tramo)
            nx = -vy
            ny = vx

            # Fuerza total normal al tramo
            Fi = q_lineal * L
            fx = Fi * nx
            fy = Fi * ny

        # Distribuir mitad a cada nodo
        nodal_forces[ni.id] += np.array([fx / 2, fy / 2])
        nodal_forces[nj.id] += np.array([fx / 2, fy / 2])


    # Paso 5: aplicar fuerzas a la estructura
    for node in nodos:
        fx, fy = nodal_forces[node.id]
        dof_x, dof_y = node.dofs
        estructura.apply_force(dof_x, fx)
        estructura.apply_force(dof_y, fy)

def apply_distributed_force_x(grupo_nodos, fuerza_total_x, estructura):
    nodos = grupo_nodos
    n = len(nodos)
    if n < 2:
        print("Se requieren al menos dos nodos para aplicar fuerza distribuida.")
        return

    # Paso 2: calcular longitud total de la línea
    longitudes = []
    total_length = 0
    for i in range(n - 1):
        dx = nodos[i + 1].x - nodos[i].x
        dy = nodos[i + 1].y - nodos[i].y
        L = np.sqrt(dx**2 + dy**2)
        longitudes.append(L)
        total_length += L

    # Paso 3: calcular carga distribuida por unidad de longitud
    q_lineal = fuerza_total_x / total_length

    # Paso 4: aplicar cargas parciales a cada nodo (2 nodos por segmento)
    nodal_forces = {node.id: np.array([0.0, 0.0]) for node in nodos}

    for i in range(n - 1):
        ni = nodos[i]
        nj = nodos[i + 1]
        xi, yi = ni.x, ni.y
        xj, yj = nj.x, nj.y

        dx = xj - xi
        dy = yj - yi
        L = longitudes[i]

        if abs(dy) < 1e-8:  # Tramo casi horizontal
            # Aplica carga puramente en X
            fx = -q_lineal * L
            fy = 0.0
        else:
            # Dirección normalizada del tramo
            vx = dx / L
            vy = dy / L

            # Vector perpendicular (hacia la izquierda respecto al tramo)
            nx = -vy
            ny = vx

            # Fuerza total normal al tramo
            Fi = q_lineal * L
            fx = Fi * nx
            fy = Fi * ny

        # Distribuir mitad a cada nodo
        nodal_forces[ni.id] += np.array([fx / 2, fy / 2])
        nodal_forces[nj.id] += np.array([fx / 2, fy / 2])

    # Paso 5: aplicar fuerzas a la estructura
    for node in nodos:
        fx, fy = nodal_forces[node.id]
        dof_x, dof_y = node.dofs
        estructura.apply_force(dof_x, fx)
        estructura.apply_force(dof_y, fy)

def apply_self_weight(elements, rho, estructure):
    P = 0
    for element in elements:
        centroid = element.get_centroid()

        area = element.area
        espesor = element.section.thickness
        peso = area * espesor * rho * 9.81 
        P += peso

        #Agrego la fuerza al elemento
        F_interna = element.apply_point_body_force(x=centroid[0], y=centroid[1],force_vector=[0, peso])

        for i in range(len(F_interna)):
            F_interna[i] =abs(F_interna[i])*-1

        node_a = element.node_list[0]
        node_b = element.node_list[1]
        node_c = element.node_list[2]

        dof_a = node_a.id * 2
        dof_b = node_b.id * 2
        dof_c = node_c.id * 2

        #Ahora agrego las cargas de peso propio

        estructure.apply_force(dof_index=dof_a, value=F_interna[1])
        estructure.apply_force(dof_index=dof_b, value=F_interna[3])
        estructure.apply_force(dof_index=dof_c, value=F_interna[5])

    print(f"Peso total de la estructura: {P:.3f} N")
    return P

def optimize_topology_iterative_n_extremes(P, grupos, elements, nodes, rho, estructure,
                                           num_iterations=5, num_elements=2,
                                           delta_t=0.2, t_min=0.2, t_max=10.0, E=0, nu=0):

    """
    Optimización topológica iterativa con propagación ultra-suavizada:
    - Aplica cambios principales a los N elementos extremos
    - Propaga ajustes suaves a través de una función Gaussiana acumulativa
    """
    import numpy as np
    import math
    from collections import defaultdict

    g = 9.81

    def gaussian_weight(level, sigma=2.0):
        return math.exp(-0.5 * (level / sigma) ** 2)

    def find_neighbors_recursive(start_indices, levels):
        neighbor_levels = defaultdict(set)
        current = set(start_indices)
        visited = set(start_indices)

        for level in range(1, levels + 1):
            next_neighbors = set()
            target_nodes = set(n for idx in current for n in elements[idx].node_list)

            for i, elem in enumerate(elements):
                if i in visited:
                    continue
                if any(n in target_nodes for n in elem.node_list):
                    neighbor_levels[level].add(i)
                    next_neighbors.add(i)

            visited.update(next_neighbors)
            current = next_neighbors

        return neighbor_levels

    def update_element_thickness(elem, delta, tag):
        t_old = elem.section.thickness
        t_new = np.clip(t_old + delta, t_min, t_max)
        elem.section = Section(t_new, E, nu)
        elem.Ke = elem.get_stiffness_matrix()
        updated_indices.add(elem.element_tag)
        #print(f"{tag} → Elem {elem.element_tag} | t: {t_old:.3f} → {t_new:.3f}")

    for it in range(num_iterations):
        print(f"\n🔁 Iteración {it+1}/{num_iterations}")
        print(f"El peso original es: {P:.5f} N")

        estructure = Solve(nodes, elements)
        apply_self_weight(elements, rho, estructure)
        #===========================================
        #AGREGAR FUERZAS DISTRIBUIDAS
        #grupo_nodos = grupos['Fuerza Y']
        #fuerza_total_y = 1200 #N NO N/m
        #apply_distributed_force_y(grupo_nodos, fuerza_total_y, estructure)
        #===========================================

        estructure.solve()

        for node in estructure.nodes:
            node.structure = estructure

        von_mises = np.array([elem.von_mises_stress(estructure.u_global) for elem in elements])
        sorted_indices = np.argsort(von_mises)

        max_indices = sorted_indices[-num_elements:]
        min_indices = sorted_indices[:num_elements]

        updated_indices = set()

        # Aplicar cambio principal
        for idx in max_indices:
            update_element_thickness(elements[idx], +delta_t, "🔺 max")

        for idx in min_indices:
            update_element_thickness(elements[idx], -delta_t, "🔻 min")

        # Propagación ultra-suavizada
        sigma = 2.0
        levels = 6  # hasta vecinos de 6º orden

        max_neighbors_by_level = find_neighbors_recursive(max_indices, levels)
        min_neighbors_by_level = find_neighbors_recursive(min_indices, levels)

        for level in range(1, levels + 1):
            weight = gaussian_weight(level, sigma) * delta_t
            for idx in max_neighbors_by_level[level]:
                if elements[idx].element_tag in updated_indices:
                    continue
                update_element_thickness(elements[idx], +weight, f"⤴ nivel {level}")
            for idx in min_neighbors_by_level[level]:
                if elements[idx].element_tag in updated_indices:
                    continue
                update_element_thickness(elements[idx], -weight, f"⤵ nivel {level}")

        # Reportar peso
        peso_total = sum(
            (el.area) * (el.section.thickness) * rho * g
            for el in elements
        )
        print(f"⚖️ Peso total aproximado: {peso_total:.5f} N")

        # Suavizar reducción de masa si el peso excede
        if peso_total > P:
            exceso = peso_total - P
            print(f"❌ Exceso de masa: {exceso:.3f} N — se reducirá espesor suavemente")

            sigma_red = 1.0
            levels_red = 6
            reduction_step = delta_t / 2
            weight_by_level = {
                level: gaussian_weight(level, sigma_red) * reduction_step
                for level in range(1, levels_red + 1)
            }

            still_exceeds = True
            temp_updated = set(updated_indices)

            max_mass_reduction_steps = 100

            all_reducible_candidates = [
                    (i, von_mises[i])
                    for i in range(len(elements))
                    if elements[i].section.thickness > t_min+0.2 and i not in temp_updated
                ]

            for reduction_iter in range(max_mass_reduction_steps):
                print(f"\n🔁 Paso de reducción #{reduction_iter + 1} — Peso actual: {peso_total:.3f} N")

                # 🔁 Recalcular esfuerzos actualizados después de cambios de espesor
                estructure = Solve(nodes, elements)
                apply_self_weight(elements, rho, estructure)
                #===========================================
                #AGREGAR FUERZAS DISTRIBUIDAS
                #grupo_nodos = grupos['Fuerza Y']
                #fuerza_total_y = 1200 #N NO N/m
                #apply_distributed_force_y(grupo_nodos, fuerza_total_y, estructure)
                #===========================================
                estructure.solve()
                von_mises = np.array([elem.von_mises_stress(estructure.u_global) for elem in elements])

                # Ordenar todos los candidatos posibles por von Mises (de menor a mayor)
                
                sorted_reduction_indices = [i for i, _ in sorted(all_reducible_candidates, key=lambda x: x[1])]

                print(f"   ➤ Candidatos a reducir: {len(sorted_reduction_indices)}")

                # Buscar `num_elements` candidatos con al menos 0.2 mm de margen
                base_indices = []
                for i in sorted_reduction_indices:
                    if elements[i].section.thickness >= t_min + 1:
                        base_indices.append(i)
                    if len(base_indices) == num_elements:
                        break

                # Si no hay suficientes con margen, completar con otros que tengan t > t_min
                if len(base_indices) < num_elements:
                    for i in sorted_reduction_indices:
                        if i in base_indices:
                            continue
                        if elements[i].section.thickness > t_min+1:
                            base_indices.append(i)
                        if len(base_indices) == num_elements:
                            break

                # Verificación final
                if not base_indices:
                    print("⚠️ No quedan elementos con espesor suficiente para seguir reduciendo masa.")
                    break

                # Aplicar reducción
                for base_idx in base_indices:
                    elem = elements[base_idx]
                    t_old = elem.section.thickness
                    t_new = max(t_old - reduction_step, t_min)
                    if t_new < t_old:
                        elem.section = Section(t_new, E, nu)
                        elem.Ke = elem.get_stiffness_matrix()
                        peso_total -= (elem.area) * ((t_old - t_new)) * rho * g
                        temp_updated.add(base_idx)

                    # Reducir vecinos suavemente
                    neighbors_by_level = find_neighbors_recursive([base_idx], levels_red)
                    for level, idxs in neighbors_by_level.items():
                        for idx in idxs:
                            if idx in temp_updated:
                                continue
                            elem_n = elements[idx]
                            if elem_n.section.thickness <= t_min:
                                continue
                            t_old_n = elem_n.section.thickness
                            t_new_n = max(t_old_n - weight_by_level[level], t_min)
                            if t_new_n < t_old_n:
                                elem_n.section = Section(t_new_n, E, nu)
                                elem_n.Ke = elem_n.get_stiffness_matrix()
                                # Recalcular peso total exacto después de todos los cambios
                                peso_total = sum(
                                    (el.area) * (el.section.thickness) * rho * g
                                    for el in elements
                                )

                                temp_updated.add(idx)

                print(f"   ➤ Nuevo peso total: {peso_total:.5f} N")

                if peso_total <= P:
                    print(f"✅ Peso ajustado suavemente en {reduction_iter + 1} pasos.")
                    still_exceeds = False
                    break

            if still_exceeds:
                print("⚠️ No fue posible ajustar completamente el peso: límite mínimo de espesor alcanzado.")

    return estructure

def compute_nodal_von_mises(elements, u_global):
        """
        Promedia los esfuerzos de Von Mises en los nodos a partir de los elementos vecinos.
        """
        nodal_vm = {}  # node.id : [list of vm from attached elements]

        for elem in elements:
            vm = elem.von_mises_stress(u_global)
            for node in elem.node_list:
                if node.id not in nodal_vm:
                    nodal_vm[node.id] = []
                nodal_vm[node.id].append(vm)

        # Promedio por nodo
        nodal_vm_avg = {node_id: np.mean(vms) for node_id, vms in nodal_vm.items()}
        return nodal_vm_avg

def main (title, lc=10, self_weight=False, distribuited_force_x = False, distribuited_force_y = False, Topologic_Optimization = False, def_scale=1, force_scale=1e-2, reaction_scale=1e-2):
    
    input_file = "PRUEBA_1_2023/geo.geo"#"VIGA_DOBLE_T/doble_t.geo"
    output_file ="PRUEBA_1_2023/msh.msh"# "VIGA_DOBLE_T/doble_t.msh"

    rho =  2500e-9#densidad, Kg/mm3 OJOOOOOO
    sigma_y_compression = 30 #MPa
    sigma_y_tension = 3#Mpa
    E = 35000#Mpa
    nu = 0.2#coeficiente de poisson

    if not (self_weight or distribuited_force_x or distribuited_force_y):
        raise ValueError("Al menos una de las opciones de carga debe ser True.")
    
    generate_mesh(input_file, output_file, lc)

    #Defino el nombre de los grupos restringidos en X e Y
    restrxy = ["Restr 3", "Restr 1", "Restr 2"]#['Restr XY 1', ....] #Si no hay hay que poner None
    restrx = None#[]
    restry = None

    grupos, mesh = make_nodes_groups(output_file, title, restrxy=restrxy, restrx=restrx, restry=restry)

    #Defino el espesor de las secciones
    thickness = {'1':3}#{'1': 115, '2': 120, '3': 120, '4': 120}, tiene que ser con numeros de Plane Surface
    sections, nodes_dict = make_sections(grupos, thickness, E, nu)
    elements, nodes = make_cst_elements(mesh, sections, nodes_dict)

    plot_all_elements(elements, title)

    estructure = Solve(nodes, elements) 

    if self_weight:
        # Aplicar peso propio a los elementos
        P = apply_self_weight(elements, rho, estructure)

    if distribuited_force_x:
        # Aplicar fuerza distribuida a los nodos
        #grupo_nodos = grupos['Fuerza X']
        #fuerza_total_x = 1200 #N NO N/m
        #apply_distributed_force_x(grupo_nodos, fuerza_total_x, estructure)
        #grupo_nodos = grupos['Fuerza X 2']
        #fuerza_total_x = 1200 #N NO N/m
        #apply_distributed_force_x(grupo_nodos, fuerza_total_x, estructure)
        pass

    if distribuited_force_y:
        # Aplicar fuerza distribuida a los nodos
        #grupo_nodos = grupos['Fuerza Y']
        #fuerza_total_y = 1200 #N NO N/m
        #apply_distributed_force_y(grupo_nodos, fuerza_total_y, estructure)
        #grupo_nodos = grupos['Fuerza Y 2']
        #fuerza_total_y = 120 #N NO N/m
        #apply_distributed_force_y(grupo_nodos, fuerza_total_y, estructure)
        pass

    start_time = time.time()  # Guarda el tiempo inicial
    estructure.solve()
    end_time = time.time()  # Guarda el tiempo después de ejecutar
    elapsed_time = end_time - start_time  # Calcula la diferencia

    # Importante: guardar los desplazamientos en cada nodo
    for node in estructure.nodes:
        node.structure = estructure  # para acceder a u_global desde cada nodo

    vm_nodal = compute_nodal_von_mises(estructure.elements, estructure.u_global)

    vm_nodal = np.array(list(vm_nodal.values()))
    
    return vm_nodal.max(), elapsed_time
    
if __name__ == "__main__":
    title = 'Convergencia' #Debe estar gentro de la carpeta GRAFICOS
    LC = [20000,10000,4000,2000, 1500, 1000, 500, 200, 100, 50, 20] #Longitud de caracteristica
    Von_mises = []
    von_mises = []
    Tiempo = []

    for lc in LC:
        print('------------------------')
        print(f'LC: {lc}')
        print('------------------------')
        resultados = main(title, lc, self_weight=True, distribuited_force_x = False, distribuited_force_y = False, Topologic_Optimization=False, def_scale=1000, force_scale=1e-2, reaction_scale=5e-3)
        vm = resultados[0]
        tiempo = resultados[1]
        Von_mises.append(vm)
        Tiempo.append(tiempo)

    
    for lc, vm in zip(LC, Von_mises):
        von_mises.append(vm)


    fig, ax1 = plt.subplots()

    # Primer eje Y: Von Mises
    ax1.set_xlabel('Characteristic Length (LC)')
    ax1.set_ylabel('Von Mises Stress (MPa)', color='tab:red')
    ax1.set_xscale('log')  # Logarithmic scale on X axis
    ax1.plot(LC, von_mises, color='tab:red', marker='o', label='Von Mises')
    ax1.tick_params(axis='y', labelcolor='tab:red')

    # Segundo eje Y: Tiempo
    ax2 = ax1.twinx()
    ax2.set_ylabel('Execution Time (s)', color='tab:blue')
    ax2.plot(LC, Tiempo, color='tab:blue', marker='s', label='Time')
    ax2.tick_params(axis='y', labelcolor='tab:blue')

    plt.title('Convergence of Von Mises Stress and Execution Time vs LC')
    fig.tight_layout()
    plt.grid(True)
    plt.savefig('GRAFICOS/convergence.png', dpi=300)
    



