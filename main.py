import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
from shapely.geometry import LineString, Point
from shapely.ops import unary_union

def Assembly(conec: list[list], C: list):
    conec = np.array(conec) - 1 # estamos convertendo conec em uma matriz numpy e fazendo decrementações em cada um dos indices para ficar mais pythonic
    nv = conec.max() + 1 #pegando valor maximo dos nos e somando 1 por conta do indice reduzido anteriormente
    nc = len(conec) #simplesmente o tanto de linhas na matriz conec
    A = np.zeros(shape=(nv,nv)) #criando matriz com valores 0 em todas as posições
    for k in range(nc):
        n1 = conec[k,0] #n1 e n2 recebem os nós de cada cano, cano=k ( de 0 até numero de canos-1 )
        n2 = conec[k,1]
        
        A[n1,n1] += C[k]
        A[n1,n2] -= C[k]
        A[n2,n1] -= C[k]
        A[n2,n2] += C[k]

    return A

def SolveNetwork(conec: list[list], C:list, natm, nB, QB):

    natm = natm - 1
    nB = nB - 1
    
    Atilde = Assembly(conec,C) #recebe a matriz A
    
    Atilde[natm, :] = 0 #todas as colunas com linha natm recebem 0
    Atilde[natm, natm] = 1 

    b = np.zeros(Atilde.shape(0)) #shape retorna uma tupla com o numero de linhas e colunas, nesse caso o 0 pega o primeiro, o numero de linhas
    b[nB]=QB
    
    pressure = np.linalg.solve(Atilde, b)
    return pressure

def PlotaRede(conec, Xno, p, q, factor_units=0.001):

    edges = conec
    coord = Xno
    nv = np.max(np.max(conec))+1
    nc = conec.shape[0]

    # Internal: get edge and midpoint coordinates
    segs = []
    mids = []
    for (i, j) in edges:
      x1, y1 = coord[i,0], coord[i,1]
      x2, y2 = coord[j,0], coord[j,1]
      segs.append(((x1, y1), (x2, y2)))
      mids.append(((x1 + x2) / 2.0, (y1 + y2) / 2.0))

    segs = np.array(segs)
    mids = np.array(mids)

    fig_size=(10, 10)
    cmap_name="coolwarm"
    node_size=50
    show_flux_labels=False
    arrow_scale=0.0005
    text_scale=0.6
    save_path=None
    show_pressure_values=False

    fig, ax = plt.subplots(figsize=fig_size)
    plt.axis('equal')

    # ---- Pressure colormap ----
    cmap = plt.get_cmap(cmap_name)
    vmin, vmax = float(p.min()), float(p.max())
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    xs, ys = [], []
    for i in range(nv):
      xs.append(coord[i,0])
      ys.append(coord[i,1])
    colors = [cmap(norm(pi)) for pi in p]
    ax.scatter(xs, ys, s=node_size, c=colors, zorder=3, edgecolors="black")

    # ---- Draw black edges and arrows ----
    #segs, mids = edge_coords()
    for idx, ((x1, y1), (x2, y2)) in enumerate(segs):
        ax.plot([x1, x2], [y1, y2], color="black", linewidth=0.5, zorder=1)

        xm, ym = mids[idx]
        dx, dy = x2 - x1, y2 - y1
        L = np.hypot(dx, dy)
        if L == 0:
          continue
        dxn, dyn = dx / L, dy / L
        nx, ny = -dyn, dxn  # normal vector

        # --- Flux arrow (black) ---
        p1, p2 = p[edges[idx,0]], p[edges[idx,1]]
        q_dir = 1 if p1 > p2 else -1

        ax.annotate(
              "",
              xy=(xm + q_dir * 0.5 * arrow_scale * dxn, ym + q_dir * 0.5 * arrow_scale * dyn),
              xytext=(xm - q_dir * 0.5 * arrow_scale * dxn, ym - q_dir * 0.5 * arrow_scale * dyn),
              arrowprops=dict(
              arrowstyle="-|>",
              color="black",
              lw=0.5,
              mutation_scale=5 * text_scale * 3,  # scales arrowhead size
              ),
              zorder=5,
              )

        # --- Flux label ---
        if show_flux_labels:
          label_offset = 0.0725*factor_units
          ax.text(
              xm + nx * label_offset,
              ym + ny * label_offset,
              f"q={q[idx]:.1e}",
              ha="center",
              va="center",
              fontsize=12 * text_scale,
              zorder=6,
              bbox=dict(facecolor='white', alpha=0.6, edgecolor='none', pad=1.0)
              )

    # ---- Node labels ----
    for node, (x, y) in enumerate(coord):
      if show_pressure_values:
          ax.text(x, y, str(node),
                  ha="center", va="center", fontsize=11 * text_scale, zorder=4)
        
          ax.text(x - 0.075*factor_units, y - 0.075*factor_units, f"p={p[node]:.1e}",
          ha="right", va="bottom", fontsize=12 * text_scale,
          color="black", zorder=5)

    # ---- Final adjustments ----
    ax.set_aspect("equal")
    ax.axis("off")

    # ---- Set limits with small margin ----
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    x_range, y_range = x_max - x_min, y_max - y_min
    ax.set_xlim(x_min - 0.5*factor_units, x_max + 0.5*factor_units)
    ax.set_ylim(y_min - 0.5*factor_units, y_max + 0.5*factor_units)

    # Optionally adapt figure size to graph geometry
    aspect_ratio = x_range / y_range if y_range != 0 else 1.0
    base_size = 8  # base figure size
    fig.set_size_inches(base_size * aspect_ratio, base_size)

    # ---- Colorbar ----
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = plt.colorbar(sm, ax=ax, label="Pressure, p", fraction=0.0225, pad=0.025)
    cbar.ax.tick_params(labelsize=10 * text_scale)
    cbar.set_label("Pressure, p", fontsize=12 * text_scale)

    if save_path:
      plt.savefig(save_path, dpi=300, bbox_inches="tight")
      plt.show()

    return fig, ax


def GeraGrafo(levels=3):
    nodes_data = []
    edges_raw = []
    node_id = 0
    
    spine_length = 6 

    # 1. Geração da Estrutura Base
    spine_nodes = []
    for i in range(spine_length):
        nodes_data.append({'x': i * 4.0, 'y': 0.0})
        spine_nodes.append(node_id)
        if i > 0: edges_raw.append((node_id - 1, node_id))
        node_id += 1
        
    def add_fractal_branches(parent_id, px, py, angle, length, depth):
        nonlocal node_id
        if depth == 0: return
        angles = [angle + np.pi/6, angle - np.pi/6]
        branch_len = length * 0.75 
        for a in angles:
            nx = px + branch_len * np.cos(a)
            ny = py + branch_len * np.sin(a)
            curr_id = node_id
            nodes_data.append({'x': nx, 'y': ny})
            edges_raw.append((parent_id, curr_id))
            node_id += 1
            add_fractal_branches(curr_id, nx, ny, a, branch_len, depth - 1)

    for s_id in spine_nodes[1:-1]:
        add_fractal_branches(s_id, nodes_data[s_id]['x'], nodes_data[s_id]['y'], np.pi/2, 3.0, levels)
        add_fractal_branches(s_id, nodes_data[s_id]['x'], nodes_data[s_id]['y'], -np.pi/2, 3.0, levels)

    # 2. Adição das Coletoras (Manifolds)
    df_temp = pd.DataFrame(nodes_data)
    y_max, y_min = df_temp['y'].max() + 1.0, df_temp['y'].min() - 1.0
    
    all_indices = [e[0] for e in edges_raw] + [e[1] for e in edges_raw]
    counts = pd.Series(all_indices).value_counts()
    leaf_ids = counts[counts == 1].index.tolist()
    leaf_ids = [idx for idx in leaf_ids if idx not in [spine_nodes[0], spine_nodes[-1]]]

    for l_id in leaf_ids:
        target_y = y_max if nodes_data[l_id]['y'] > 0 else y_min
        new_id = node_id
        nodes_data.append({'x': nodes_data[l_id]['x'], 'y': target_y})
        edges_raw.append((l_id, new_id))
        node_id += 1

    # 3. Processamento Geométrico de Interseções
    lines = [LineString([(nodes_data[e[0]]['x'], nodes_data[e[0]]['y']), 
                         (nodes_data[e[1]]['x'], nodes_data[e[1]]['y'])]) for e in edges_raw]
    
    # Adicionar linhas das coletoras para o merge
    df_nodes_final = pd.DataFrame(nodes_data)
    for y_lim in [y_max, y_min]:
        pts = df_nodes_final[df_nodes_final['y'] == y_lim].sort_values('x')
        if len(pts) > 1:
            lines.append(LineString(pts[['x', 'y']].values))

    # Quebra todas as linhas nas interseções
    merged_graph = unary_union(lines)

    # 4. Conversão para Arrays NumPy (Mapeamento de IDs)
    final_nodes_map = {}
    final_nodes_list = []
    final_edges_list = []
    
    def get_node_id(pt):
        # Arredondamento para evitar erros de precisão de ponto flutuante
        coords = (round(pt[0], 6), round(pt[1], 6))
        if coords not in final_nodes_map:
            final_nodes_map[coords] = len(final_nodes_list)
            final_nodes_list.append([pt[0], pt[1]])
        return final_nodes_map[coords]

    # Itera sobre cada segmento gerado pelo unary_union
    segments = merged_graph.geoms if hasattr(merged_graph, 'geoms') else [merged_graph]
    for seg in segments:
        id_start = get_node_id(seg.coords[0])
        id_end = get_node_id(seg.coords[-1])
        final_edges_list.append([id_start, id_end])

    return np.array(final_nodes_list), np.array(final_edges_list)