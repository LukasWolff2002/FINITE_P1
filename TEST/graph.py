import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import os
import re

def clean_filename(text):
    """Convierte un string en nombre de archivo válido eliminando LaTeX y símbolos especiales."""
    return re.sub(r'[^\w\-_. ]', '', text.replace('$', '').replace('\\', '').replace('{', '').replace('}', ''))

def plot_all_elements(elements, title, show_ids=True):
    all_x = []
    all_y = []

    # Recopilar coordenadas antes de crear la figura
    for elem in elements:
        coords = elem.get_xy_matrix()
        coords = np.vstack([coords, coords[0]])
        all_x.extend(coords[:, 0])
        all_y.extend(coords[:, 1])

    # Definir márgenes y límites
    x_min, x_max = min(all_x), max(all_x)
    y_min, y_max = min(all_y), max(all_y)
    x_margin = (x_max - x_min) * 0.05
    y_margin = (y_max - y_min) * 0.05

    x_range = (x_max - x_min) + 2 * x_margin
    y_range = (y_max - y_min) + 2 * y_margin

    # Fijar ancho (en pulgadas), y escalar alto proporcionalmente
    fixed_width = 8  # Puedes cambiarlo
    aspect_ratio = y_range / x_range
    height = fixed_width * aspect_ratio

    fig, ax = plt.subplots(figsize=(fixed_width, height))

    # Dibujar elementos
    for elem in elements:
        coords = elem.get_xy_matrix()
        coords = np.vstack([coords, coords[0]])
        ax.plot(coords[:, 0], coords[:, 1], 'k-', linewidth=1)

        if show_ids:
            for nodo, (x, y) in zip(elem.node_list, coords[:-1]):
                if np.all(nodo.restrain == [1, 1]):
                    pass
                    # ax.text(x, y, f'N{nodo.id}', color='blue', fontsize=6, ha='center', va='center')

    # Ajustar límites y aspecto
    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    ax.set_ylim(y_min - y_margin, y_max + y_margin)
    ax.set_aspect('equal', adjustable='box')

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("All CST elements")

    ax.grid(True)

    # Guardar imagen recortando al contenido real
    plt.savefig(f"GRAFICOS/{title}_elementos.png", dpi=300, bbox_inches='tight')
    plt.close()


def plot_applied_forces(nodes, elements, title, f_vector, scale=1e-2):
    all_x = []
    all_y = []

    # Recolectar coordenadas de todos los elementos
    for elem in elements:
        coords = elem.get_xy_matrix()
        all_x.extend(coords[:, 0])
        all_y.extend(coords[:, 1])

    # También considerar nodos, por si hay fuerzas fuera del dominio de elementos
    for node in nodes:
        all_x.append(node.x)
        all_y.append(node.y)

    # Calcular límites con márgenes
    x_min, x_max = min(all_x), max(all_x)
    y_min, y_max = min(all_y), max(all_y)
    x_margin = (x_max - x_min) * 0.05
    y_margin = (y_max - y_min) * 0.05

    x_range = (x_max - x_min) + 2 * x_margin
    y_range = (y_max - y_min) + 2 * y_margin

    # Escalar altura según ancho fijo
    fixed_width = 8  # pulgadas
    aspect_ratio = y_range / x_range
    height = fixed_width * aspect_ratio

    fig, ax = plt.subplots(figsize=(fixed_width, height))
    ax.set_title("Forces applied on nodes")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    ax.set_ylim(y_min - y_margin, y_max + y_margin)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True)

    # Dibujar elementos
    for elem in elements:
        coords = elem.get_xy_matrix()
        coords = np.vstack([coords, coords[0]])
        ax.plot(coords[:, 0], coords[:, 1], color='lightgray', linewidth=1, zorder=1)

    # Dibujar fuerzas
    for node in nodes:
        dof_x, dof_y = node.dofs
        fx = f_vector[dof_x][0] if dof_x < len(f_vector) else 0.0
        fy = f_vector[dof_y][0] if dof_y < len(f_vector) else 0.0

        if not np.isclose(fx, 0.0) or not np.isclose(fy, 0.0):
            ax.quiver(node.x, node.y, fx, fy,
                      angles='xy', scale_units='xy',
                      scale=1 / scale, color='red', width=0.005,
                      zorder=5)

    # Guardar sin márgenes blancos
    plt.savefig(f"GRAFICOS/{title}_fuerzas.png", dpi=300, bbox_inches='tight')
    plt.close()

import matplotlib.tri as mtri

def plot_deformed_structure(elements, title, scale=1.0, show_ids=False):
    # Listas para las coordenadas deformadas y sus desplazamientos
    points = []
    displacements = []
    triangles = []
    point_index = {}
    idx_counter = 0

    for elem in elements:
        coords = elem.get_xy_matrix()
        u = np.zeros_like(coords)

        # Obtener desplazamientos globales
        for i, nodo in enumerate(elem.node_list):
            ux = nodo.structure.u_global[nodo.dofs[0], 0]
            uy = nodo.structure.u_global[nodo.dofs[1], 0]
            u[i] = [ux, uy]

        # Coordenadas deformadas
        def_coords = coords + scale * u

        tri = []
        for i in range(3):
            key = tuple(def_coords[i])
            if key not in point_index:
                point_index[key] = idx_counter
                points.append(key)
                # Magnitud del desplazamiento en ese nodo
                mag = np.linalg.norm(u[i])
                displacements.append(mag)
                idx_counter += 1
            tri.append(point_index[key])

        triangles.append(tri)

    points = np.array(points)
    displacements = np.array(displacements)
    triangles = np.array(triangles)

    # Crear triangulación y gráfico
    fig, ax = plt.subplots(figsize=(8, 6))
    triang = mtri.Triangulation(points[:, 0], points[:, 1], triangles)

    tpc = ax.tripcolor(triang, displacements/10, shading='gouraud', cmap='viridis')
    cb = plt.colorbar(tpc, ax=ax)
    cb.set_label('Displacement magnitude cm')

    ax.set_aspect('equal')
    ax.set_title(f"Estructura deformada (×{scale}) - Mapa de calor")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True)

    plt.savefig(f"GRAFICOS/{title}_deformada_colormap.png", dpi=300, bbox_inches='tight')
    plt.close()

def plot_deformed_with_reactions(title, elements, reactions, scale=1.0, reaction_scale=1e-3, show_ids=False):
    import matplotlib.tri as mtri
    all_x = []
    all_y = []

    for elem in elements:
        coords = elem.get_xy_matrix()
        all_x.extend(coords[:, 0])
        all_y.extend(coords[:, 1])

    x_min, x_max = min(all_x), max(all_x)
    y_min, y_max = min(all_y), max(all_y)
    x_margin = (x_max - x_min) * 0.05
    y_margin = (y_max - y_min) * 0.05

    x_range = (x_max - x_min) + 2 * x_margin
    y_range = (y_max - y_min) + 2 * y_margin

    fixed_width = 16
    aspect_ratio = y_range / x_range
    height = fixed_width * aspect_ratio / 2

    fig, axs = plt.subplots(1, 2, figsize=(fixed_width, height))
    ax1, ax2 = axs

    # === GRÁFICO 1: Reacciones solamente ===
    for elem in elements:
        coords = elem.get_xy_matrix()
        coords_closed = np.vstack([coords, coords[0]])
        ax1.plot(coords_closed[:, 0], coords_closed[:, 1], 'lightgray', linewidth=0.5)

        for i, nodo in enumerate(elem.node_list):
            rx = reactions[nodo.dofs[0]][0] if nodo.restrain[0] == 1 else 0.0
            ry = reactions[nodo.dofs[1]][0] if nodo.restrain[1] == 1 else 0.0

            if rx != 0.0 or ry != 0.0:
                ax1.quiver(nodo.x, nodo.y, rx, ry, angles='xy', scale_units='xy',
                           scale=1/reaction_scale, color='red', width=0.005)
            if show_ids:
                ax1.text(nodo.x, nodo.y, f'N{nodo.id}', fontsize=6)

    ax1.set_title("Nodal reactions", fontsize=14)
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.grid(True)
    ax1.set_xlim(x_min - x_margin, x_max + x_margin)
    ax1.set_ylim(y_min - y_margin, y_max + y_margin)
    ax1.set_aspect('equal', adjustable='box')

    # === GRÁFICO 2: Heatmap de magnitud de reacciones ===
    nodal_values = {}
    for elem in elements:
        for nodo in elem.node_list:
            rx = reactions[nodo.dofs[0]][0] if nodo.restrain[0] == 1 else 0.0
            ry = reactions[nodo.dofs[1]][0] if nodo.restrain[1] == 1 else 0.0
            mag = np.sqrt(rx**2 + ry**2)
            nodal_values[nodo.id] = mag

    all_nodes = {n.id: n for e in elements for n in e.node_list}
    node_list = list(all_nodes.values())
    node_id_map = {node.id: i for i, node in enumerate(node_list)}
    points = np.array([[node.x, node.y] for node in node_list])
    triangles = []
    values = []

    for elem in elements:
        ids = [n.id for n in elem.node_list]
        val = np.mean([nodal_values.get(nid, 0.0) for nid in ids])
        if val > 0:
            tri = [node_id_map[nid] for nid in ids]
            triangles.append(tri)
            values.append(val)

    triang = mtri.Triangulation(points[:, 0], points[:, 1], triangles)

    for elem in elements:
        coords = elem.get_xy_matrix()
        coords = np.vstack([coords, coords[0]])
        ax2.plot(coords[:, 0], coords[:, 1], color='lightgray', linewidth=0.5)

    if values:
        tpc = ax2.tripcolor(triang, facecolors=values, edgecolors='k', cmap='plasma', zorder=2)
        cbar = fig.colorbar(tpc, ax=ax2)
        cbar.set_label("Reaction force [N]", fontsize=14)

    ax2.set_title("Heatmap of reactions per element", fontsize=14)
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    ax2.grid(True)
    ax2.set_xlim(x_min - x_margin, x_max + x_margin)
    ax2.set_ylim(y_min - y_margin, y_max + y_margin)
    ax2.set_aspect('equal', adjustable='box')

    fig.savefig(f"GRAFICOS/{title}_deformada_reacciones.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Obtener nodos únicos una sola vez
    all_nodes = {n.id: n for e in elements for n in e.node_list}.values()

    total_rx = 0.0
    total_ry = 0.0

    for node in all_nodes:
        rx = reactions[node.dofs[0], 0] if node.restrain[0] == 1 else 0.0
        ry = reactions[node.dofs[1], 0] if node.restrain[1] == 1 else 0.0

        total_rx += rx
        total_ry += ry

    total_magnitude = np.sqrt(total_rx**2 + total_ry**2)
    print(f"📌 Sumatoria total de reacciones: Rx = {total_rx:.3f} N, Ry = {total_ry:.3f} N, ||R|| = {total_magnitude:.3f} N")


def plot_von_mises_field(nodes, elements, vm_nodal_dict, title, cmap='plasma'):
    node_id_to_index = {}
    xs, ys, vms = [], [], []

    for i, node in enumerate(nodes):
        node_id_to_index[node.id] = i
        xs.append(node.x)
        ys.append(node.y)
        vms.append(vm_nodal_dict.get(node.id, 0.0))

    triangles = []
    for elem in elements:
        triangle = [node_id_to_index[n.id] for n in elem.node_list]
        triangles.append(triangle)

    triang = mtri.Triangulation(xs, ys, triangles)

    # Calcular límites para adaptar proporción del gráfico
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    x_margin = (x_max - x_min) * 0.05
    y_margin = (y_max - y_min) * 0.05

    x_range = (x_max - x_min) + 2 * x_margin
    y_range = (y_max - y_min) + 2 * y_margin

    fixed_width = 8  # pulgadas
    aspect_ratio = y_range / x_range
    height = fixed_width * aspect_ratio

    fig, ax = plt.subplots(figsize=(fixed_width, height))

    tcf = ax.tricontourf(triang, vms, levels=20, cmap=cmap)
    #ax.triplot(triang, color='gray', linewidth=0.5)
    cbar = fig.colorbar(tcf, ax=ax)
    cbar.set_label("Von Mises Stress (MPa)")

    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    ax.set_ylim(y_min - y_margin, y_max + y_margin)
    ax.set_aspect('equal', adjustable='box')
    ax.set_title("Von Mises stress field over elements")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True)

    # Guardar imagen recortada y centrada
    fig.savefig(f"GRAFICOS/{title}_von_mises.png", dpi=300, bbox_inches='tight')
    plt.close()

from matplotlib.colors import BoundaryNorm

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import LogNorm

def plot_von_mises_per_element(nodes, elements, vm_nodal_dict, title, cmap='plasma'):
    node_id_to_index = {}
    xs, ys = [], []

    for i, node in enumerate(nodes):
        node_id_to_index[node.id] = i
        xs.append(node.x)
        ys.append(node.y)

    triangles = []
    element_colors = []

    for elem in elements:
        triangle = [node_id_to_index[n.id] for n in elem.node_list]
        triangles.append(triangle)

        # Usar promedio de los nodos para el valor del elemento
        vms_elem = np.mean([vm_nodal_dict.get(n.id, 0.0) for n in elem.node_list])

        # ⚡ Si el valor es 0 o negativo, reemplazar por pequeño valor positivo
        if vms_elem <= 0:
            vms_elem = 1e-6

        element_colors.append(vms_elem)

    triang = mtri.Triangulation(xs, ys, triangles)

    # --- Escalado proporcional ---
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    x_margin = (x_max - x_min) * 0.05
    y_margin = (y_max - y_min) * 0.05

    x_range = (x_max - x_min) + 2 * x_margin
    y_range = (y_max - y_min) + 2 * y_margin

    fixed_width = 8
    aspect_ratio = y_range / x_range
    height = fixed_width * aspect_ratio

    fig, ax = plt.subplots(figsize=(fixed_width, height))

    # --- Escala logarítmica ---
    vmin = min(element_colors)
    vmax = max(element_colors)

    if vmin <= 0:
        vmin = min([v for v in element_colors if v > 0], default=1e-6)

    norm = LogNorm(vmin=vmin, vmax=vmax)

    # 🎯 Graficar con shading='flat'
    tpc = ax.tripcolor(
        triang,
        facecolors=element_colors,
        cmap=cmap,
        norm=norm,
        edgecolors='k',
        shading='flat'
    )

    cbar = fig.colorbar(tpc, ax=ax)
    cbar.set_label("Von Mises Stress (MPa)")

    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    ax.set_ylim(y_min - y_margin, y_max + y_margin)
    ax.set_aspect('equal', adjustable='box')
    ax.set_title("Von Mises stress per element")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True)

    os.makedirs("GRAFICOS", exist_ok=True)
    fig.savefig(f"GRAFICOS/{title}_von_mises_per_element.png", dpi=300, bbox_inches='tight')
    plt.close()


def compute_nodal_principal_fields(nodes, elements, u_global):
    """
    Calcula los esfuerzos y deformaciones principales σ1, σ2, ε1, ε2 por nodo.
    Devuelve 4 diccionarios: sigma1, sigma2, eps1, eps2
    """
    sigma1_map = {node.id: [] for node in nodes}
    sigma2_map = {node.id: [] for node in nodes}
    eps1_map = {node.id: [] for node in nodes}
    eps2_map = {node.id: [] for node in nodes}

    for elem in elements:
        σ = elem.get_stress(u_global)
        ε = elem.get_strain(u_global)

        σx, σy, τxy = σ
        εx, εy, γxy = ε

        # Esfuerzos principales
        σ_avg = 0.5 * (σx + σy)
        Rσ = np.sqrt(((σx - σy) / 2)**2 + τxy**2)
        σ1, σ2 = σ_avg + Rσ, σ_avg - Rσ

        # Deformaciones principales
        ε_avg = 0.5 * (εx + εy)
        Rε = np.sqrt(((εx - εy) / 2)**2 + (γxy / 2)**2)
        ε1, ε2 = ε_avg + Rε, ε_avg - Rε

        for node in elem.node_list:
            sigma1_map[node.id].append(σ1)
            sigma2_map[node.id].append(σ2)
            eps1_map[node.id].append(ε1)
            eps2_map[node.id].append(ε2)

    # Promediar por nodo
    sigma1 = {nid: np.mean(vals) for nid, vals in sigma1_map.items()}
    sigma2 = {nid: np.mean(vals) for nid, vals in sigma2_map.items()}
    eps1 = {nid: np.mean(vals) for nid, vals in eps1_map.items()}
    eps2 = {nid: np.mean(vals) for nid, vals in eps2_map.items()}

    return sigma1, sigma2, eps1, eps2

def plot_all_scalar_fields_separately(nodes, elements, nodal_fields, title_prefix):
    """
    Genera 6 gráficos individuales para: σxx, σyy, τxy, εxx, εyy, γxy
    """
    σxx, σyy, τxy, εxx, εyy, γxy = {}, {}, {}, {}, {}, {}

    # Llenar los diccionarios con los valores de las tensiones y deformaciones
    for node_id, (sigma, epsilon) in nodal_fields.items():
        σxx[node_id] = sigma[0]
        σyy[node_id] = sigma[1]
        τxy[node_id] = sigma[2]
        εxx[node_id] = epsilon[0]
        εyy[node_id] = epsilon[1]
        γxy[node_id] = epsilon[2]

    # Función para imprimir los valores máximos y mínimos
    def print_min_max(field, name):
        max_value = max(field.values())
        min_value = min(field.values())
        print(f"{name} - Max: {max_value:.2f}, Min: {min_value:.2f}")

    # Imprimir los valores máximos y mínimos de cada campo
    print_min_max(σxx, r"$\sigma_{xx}$ (MPa)")
    print_min_max(σyy, r"$\sigma_{yy}$ (MPa)")
    print_min_max(τxy, r"$\tau_{xy}$ (MPa)")
    print_min_max(εxx, r"$\epsilon_{xx}$")
    print_min_max(εyy, r"$\epsilon_{yy}$")
    print_min_max(γxy, r"$\gamma_{xy}$")

    plot_scalar_field(nodes, elements, σxx, r"$\sigma_{xx}$ (MPa)", f"{title_prefix} - sigma_xx")
    plot_scalar_field(nodes, elements, σyy, r"$\sigma_{yy}$ (MPa)", f"{title_prefix} - sigma_yy")
    plot_scalar_field(nodes, elements, τxy, r"$\tau_{xy}$ (MPa)", f"{title_prefix} - tau_xy")
    plot_scalar_field(nodes, elements, εxx, r"$\varepsilon_{xx}$", f"{title_prefix} - epsilon_xx")
    plot_scalar_field(nodes, elements, εyy, r"$\varepsilon_{yy}$", f"{title_prefix} - epsilon_yy")
    plot_scalar_field(nodes, elements, γxy, r"$\gamma_{xy}$", f"{title_prefix} - gamma_xy")

    plot_scalar_field_per_element(nodes, elements, σxx, r"$\sigma_{xx}$ (MPa)", f"{title_prefix} - sigma_xx_per_element")
    plot_scalar_field_per_element(nodes, elements, σyy, r"$\sigma_{yy}$ (MPa)", f"{title_prefix} - sigma_yy_per_element")
    plot_scalar_field_per_element(nodes, elements, τxy, r"$\tau_{xy}$ (MPa)", f"{title_prefix} - tau_xy_per_element")
    plot_scalar_field_per_element(nodes, elements, εxx, r"$\varepsilon_{xx}$", f"{title_prefix} - epsilon_xx_per_element")
    plot_scalar_field_per_element(nodes, elements, εyy, r"$\varepsilon_{yy}$", f"{title_prefix} - epsilon_yy_per_element")
    plot_scalar_field_per_element(nodes, elements, γxy, r"$\gamma_{xy}$", f"{title_prefix} - gamma_xy_per_element")

def plot_scalar_field(nodes, elements, nodal_values, field_title, filename_prefix, cmap='plasma'):
    """
    Grafica un campo escalar interpolado sobre la malla de elementos.
    """
    node_id_to_index = {node.id: i for i, node in enumerate(nodes)}
    xs = [node.x for node in nodes]
    ys = [node.y for node in nodes]

    # Triángulos de elementos
    triangles = [
        [node_id_to_index[n.id] for n in elem.node_list]
        for elem in elements
    ]

    # Valores del campo escalar por nodo
    values = np.zeros(len(nodes))
    for node in nodes:
        values[node_id_to_index[node.id]] = nodal_values[node.id]

    # Márgenes para recorte automático
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    x_margin = (x_max - x_min) * 0.05
    y_margin = (y_max - y_min) * 0.05

    x_range = (x_max - x_min) + 2 * x_margin
    y_range = (y_max - y_min) + 2 * y_margin

    fixed_width = 8  # en pulgadas
    aspect_ratio = y_range / x_range
    height = fixed_width * aspect_ratio

    fig, ax = plt.subplots(figsize=(fixed_width, height))
    triang = mtri.Triangulation(xs, ys, triangles)

    tcf = ax.tricontourf(triang, values, levels=20, cmap=cmap)
    cbar = fig.colorbar(tcf, ax=ax)
    cbar.set_label(field_title)

    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    ax.set_ylim(y_min - y_margin, y_max + y_margin)
    ax.set_aspect('equal', adjustable='box')
    ax.set_title(field_title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True)

    # Guardar imagen recortada
    os.makedirs("GRAFICOS", exist_ok=True)
    filename = f"GRAFICOS/{filename_prefix}.png"
    fig.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize

def plot_scalar_field_per_element(nodes, elements, nodal_values, field_title, filename_prefix, cmap='plasma'):
    node_id_to_index = {node.id: i for i, node in enumerate(nodes)}
    xs = [node.x for node in nodes]
    ys = [node.y for node in nodes]

    triangles = []
    element_values = []

    for elem in elements:
        try:
            triangle = [node_id_to_index[n.id] for n in elem.node_list]
            value = np.mean([nodal_values[n.id] for n in elem.node_list])

            #if value <= 0:
            #    value = 1e-6

            triangles.append(triangle)
            element_values.append(value)
        except KeyError:
            # Si un nodo no tiene valor, se ignora el triángulo
            continue

    # Conversión a numpy
    triangles = np.array(triangles)
    element_values = np.array(element_values)

    triang = mtri.Triangulation(xs, ys, triangles)

    # Dimensiones del gráfico
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    x_margin = (x_max - x_min) * 0.05
    y_margin = (y_max - y_min) * 0.05
    aspect_ratio = (y_max - y_min + 2 * y_margin) / (x_max - x_min + 2 * x_margin)
    height = 8 * aspect_ratio

    fig, ax = plt.subplots(figsize=(8, height))

    vmin = np.min(element_values)
    vmax = np.max(element_values)

    norm = Normalize(vmin=vmin, vmax=vmax)

    tpc = ax.tripcolor(
        triang,
        facecolors=element_values,
        cmap=cmap,
        norm=norm,
        edgecolors='k',
        shading='flat'
    )

    cbar = fig.colorbar(tpc, ax=ax)
    cbar.set_label(field_title)

    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    ax.set_ylim(y_min - y_margin, y_max + y_margin)
    ax.set_aspect('equal')
    ax.set_title(field_title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True)

    os.makedirs("GRAFICOS", exist_ok=True)
    fig.savefig(f"GRAFICOS/{filename_prefix}.png", dpi=300, bbox_inches='tight')
    plt.close()

from matplotlib.colors import LogNorm
import matplotlib.tri as mtri

def plot_safety_field_per_element_log(nodes, elements, nodal_values, field_title, filename_prefix, cmap='viridis_r'):
    """
    Grafica un campo de factores de seguridad (FS) por elemento usando escala logarítmica.
    """

    node_id_to_index = {node.id: i for i, node in enumerate(nodes)}
    xs = [node.x for node in nodes]
    ys = [node.y for node in nodes]

    triangles = []
    element_values = []

    for elem in elements:
        try:
            triangle = [node_id_to_index[n.id] for n in elem.node_list]
            value = np.mean([nodal_values[n.id] for n in elem.node_list])

            # Seguridad para LogNorm: reemplazar valores pequeños o negativos
            if value <= 0:
                value = 1e-3  # Muy pequeño pero válido para escala log
            triangles.append(triangle)
            element_values.append(value)
        except KeyError:
            continue

    triangles = np.array(triangles)
    element_values = np.array(element_values)

    triang = mtri.Triangulation(xs, ys, triangles)

    # Definir límites del gráfico
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    x_margin = (x_max - x_min) * 0.05
    y_margin = (y_max - y_min) * 0.05
    aspect_ratio = (y_max - y_min + 2 * y_margin) / (x_max - x_min + 2 * x_margin)
    height = 8 * aspect_ratio

    fig, ax = plt.subplots(figsize=(8, height))

    vmin = max(np.min(element_values), 1e-3)
    vmax = np.max(element_values)
    norm = LogNorm(vmin=vmin, vmax=vmax)

    tpc = ax.tripcolor(
        triang,
        facecolors=element_values,
        cmap=cmap,
        norm=norm,
        edgecolors='k',
        shading='flat'
    )

    cbar = fig.colorbar(tpc, ax=ax)
    cbar.set_label(field_title)

    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    ax.set_ylim(y_min - y_margin, y_max + y_margin)
    ax.set_aspect('equal')
    ax.set_title(field_title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True)

    os.makedirs("GRAFICOS", exist_ok=True)
    fig.savefig(f"GRAFICOS/{filename_prefix}.png", dpi=300, bbox_inches='tight')
    plt.close()


def plot_safety_factors(nodes, elements, u_global, title_prefix, sigma_y_tension, sigma_y_compression):
    """
    Genera dos mapas de factor de seguridad:
    - Tracción: FS = sigma_y_tension / abs(sigma_1)
    - Compresión: FS = sigma_y_compression / abs(sigma_2)
    También imprime y muestra el mínimo FS en cada caso.
    """
    # Obtener tensiones principales por nodo
    sigma1, sigma2, _, _ = compute_nodal_principal_fields(nodes, elements, u_global)

    # Calcular FS por nodo
    fs_tension = {nid: sigma_y_tension / max(abs(val), 1e-6) for nid, val in sigma1.items()}
    fs_compression = {nid: sigma_y_compression / max(abs(val), 1e-6) for nid, val in sigma2.items()}

    # Calcular mínimos
    min_fs_tension = min(fs_tension.values())
    min_fs_compression = min(fs_compression.values())

    print(f"🔧 Mínimo FS Tracción: {min_fs_tension:.3f}")
    print(f"🔧 Mínimo FS Compresión: {min_fs_compression:.3f}")

    # Graficar por elemento
    plot_safety_field_per_element_log(
        nodes, elements, fs_tension,
        field_title=f"FS Tracción (min = {min_fs_tension:.2f})",
        filename_prefix=f"{title_prefix}_FS_Tension_Log"
    )

    plot_safety_field_per_element_log(
        nodes, elements, fs_compression,
        field_title=f"FS Compresión (min = {min_fs_compression:.2f})",
        filename_prefix=f"{title_prefix}_FS_Compression_Log"
    )


def plot_principal_fields(nodes, elements, u_global, title_prefix, sigma_y_tension, sigma_y_compression):
        """
        Genera 4 gráficos: σ1, σ2, ε1, ε2.
        """
        sigma1, sigma2, eps1, eps2 = compute_nodal_principal_fields(nodes, elements, u_global)

        plot_scalar_field(nodes, elements, sigma1, r"$\sigma_1$ (MPa)", f"{title_prefix} - sigma_1")
        plot_scalar_field(nodes, elements, sigma2, r"$\sigma_2$ (MPa)", f"{title_prefix} - sigma_2")
        plot_scalar_field(nodes, elements, eps1, r"$\varepsilon_1$", f"{title_prefix} - epsilon_1")
        plot_scalar_field(nodes, elements, eps2, r"$\varepsilon_2$", f"{title_prefix} - epsilon_2")

        plot_scalar_field_per_element(nodes, elements, sigma1, r"$\sigma_1$ (MPa)", f"{title_prefix} - sigma_1_per_element")
        plot_scalar_field_per_element(nodes, elements, sigma2, r"$\sigma_2$ (MPa)", f"{title_prefix} - sigma_2_per_element")
        plot_scalar_field_per_element(nodes, elements, eps1, r"$\varepsilon_1$", f"{title_prefix} - epsilon_1_per_element")
        plot_scalar_field_per_element(nodes, elements, eps2, r"$\varepsilon_2$", f"{title_prefix} - epsilon_2_per_element")
        plot_safety_factors(nodes, elements, u_global,title_prefix, sigma_y_tension, sigma_y_compression)

import matplotlib.tri as mtri
import matplotlib.pyplot as plt
import numpy as np

def plot_elements_by_thickness(elements, title, cmap='viridis'):
    """
    Dibuja y guarda una figura coloreando los elementos según el espesor de su sección,
    ajustando automáticamente el alto para mantener proporciones.
    """
    node_ids = {}
    xs, ys = [], []
    triangles = []
    thicknesses = []

    counter = 0
    for elem in elements:
        triangle = []
        for node in elem.node_list:
            if node.id not in node_ids:
                node_ids[node.id] = counter
                xs.append(node.x)
                ys.append(node.y)
                counter += 1
            triangle.append(node_ids[node.id])
        triangles.append(triangle)
        thicknesses.append(elem.section.thickness)

    triang = mtri.Triangulation(xs, ys, triangles)

    # Cálculo de límites y dimensiones
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    x_margin = (x_max - x_min) * 0.05
    y_margin = (y_max - y_min) * 0.05

    x_range = (x_max - x_min) + 2 * x_margin
    y_range = (y_max - y_min) + 2 * y_margin

    fixed_width = 8
    aspect_ratio = y_range / x_range
    height = fixed_width * aspect_ratio

    fig, ax = plt.subplots(figsize=(fixed_width, height))
    tpc = ax.tripcolor(triang, facecolors=thicknesses, edgecolors='k', cmap=cmap)
    cbar = plt.colorbar(tpc, ax=ax, label="Thickness (mm)")
    ax.set_title("Thickness Distribution per Element")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    ax.set_ylim(y_min - y_margin, y_max + y_margin)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True)

    plt.tight_layout()
    plt.savefig(f"GRAFICOS/{title}_espesores_post_topo.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Gráfico guardado como GRAFICOS/{title}_espesores.png")

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

def compute_nodal_stress_strain(nodes, elements, u_global):

    """
    Calcula σxx, σyy, σxy, εxx, εyy, εxy por nodo (promedio de elementos conectados).
    Retorna: diccionario {node.id: (σ_vec, ε_vec)}
    """
    stress_map = {node.id: [] for node in nodes}
    strain_map = {node.id: [] for node in nodes}

    for elem in elements:
        stress = elem.get_stress(u_global)
        strain = elem.get_strain(u_global)
        for node in elem.node_list:
            stress_map[node.id].append(stress)
            strain_map[node.id].append(strain)

    result = {}
    for node in nodes:
        σ_avg = np.mean(stress_map[node.id], axis=0) if stress_map[node.id] else np.zeros(3)
        ε_avg = np.mean(strain_map[node.id], axis=0) if strain_map[node.id] else np.zeros(3)
        result[node.id] = (σ_avg, ε_avg)

    return result  # node.id: (σ_vec, ε_vec)

def plot_results (estructure, elements, title, def_scale=1, force_scale=1e-2, reaction_scale=1e-2, sigma_y_tension=0, sigma_y_compression=0):
    f = estructure.f_original if hasattr(estructure, 'f_original') else estructure.f_global

    plot_applied_forces(estructure.nodes, elements,title, f, scale=force_scale)

    desplazamientos = estructure.solve()  # Resuelve la estructura y obtiene los desplazamientos

    # Inicializamos la variable para almacenar el máximo desplazamiento vertical
    max_desplazamiento_vertical = 0
    min_desplazamiento_vertical = 0

    # Iteramos sobre los nodos de la estructura
    for node in estructure.nodes:
        # El índice para el desplazamiento en y (vertical) será el índice impar de los DOF
        # Asumimos que el desplazamiento en y está en el índice impar correspondiente al nodo
        desplazamiento_vertical = desplazamientos[node.dofs[1]]  # DOF[1] es el desplazamiento en y (vertical)

        # Comprobamos si este es el mayor desplazamiento vertical encontrado
        if desplazamiento_vertical < max_desplazamiento_vertical:
            max_desplazamiento_vertical = desplazamiento_vertical

        if desplazamiento_vertical > min_desplazamiento_vertical:
            min_desplazamiento_vertical = desplazamiento_vertical

        # Si deseas imprimir el desplazamiento de un nodo específico, puedes hacerlo aquí:
        if node.x == 2800 and node.y == 1600:
            print(f"Desplazamiento en nodo {node.id}: x = {node.x}, y = {node.y}")
            print(f"Desplazamiento x: {desplazamientos[node.dofs[0]][0]}, y: {desplazamientos[node.dofs[1]][0]} mm")

    # Imprimir el desplazamiento máximo vertical
    print(f"El desplazamiento máximo vertical es {max_desplazamiento_vertical} mm")
    print(f"El desplazamiento mínimo vertical es {min_desplazamiento_vertical} mm")



    # Importante: guardar los desplazamientos en cada nodo
    for node in estructure.nodes:
        node.structure = estructure  # para acceder a u_global desde cada nodo

    # Luego graficar
    plot_deformed_structure(estructure.elements, title, scale=def_scale, show_ids=False)

    reacciones = estructure.compute_reactions()

    plot_deformed_with_reactions(title, estructure.elements, reacciones, scale=def_scale, reaction_scale=reaction_scale, show_ids=False)

    vm_nodal = compute_nodal_von_mises(estructure.elements, estructure.u_global)
    plot_von_mises_field(estructure.nodes, estructure.elements, vm_nodal, title)
    plot_von_mises_per_element(estructure.nodes, estructure.elements, vm_nodal, title)
    
    nodal_fields = compute_nodal_stress_strain(estructure.nodes, estructure.elements, estructure.u_global)   
    
    plot_all_scalar_fields_separately(estructure.nodes, estructure.elements, nodal_fields, title)
    plot_principal_fields(estructure.nodes, estructure.elements, estructure.u_global, title_prefix=title, sigma_y_tension=sigma_y_tension, sigma_y_compression=sigma_y_compression)