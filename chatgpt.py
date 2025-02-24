import itertools
import networkx as nx
import matplotlib.pyplot as plt

#----------------------
# 1. 定義 C_matrix, M_matrix, Pa_names, Pa_labels
#----------------------
C_matrix = {
    'S1': [1, 0, 0, 0, 1, 0],  # S1 -> Pa1, Pa5
    'S2': [1, 0, 0, 0, 0, 1],  # S2 -> Pa1, Pa6
    'S3': [1, 1, 0, 0, 0, 0],  # S3 -> Pa1, Pa2
    'S4': [1, 0, 1, 1, 0, 0],  # S4 -> Pa1, Pa3, Pa4
    'S5': [1, 0, 1, 1, 0, 0],  # S5 -> Pa1, Pa3, Pa4
    'S6': [1, 1, 0, 1, 0, 0]   # S6 -> Pa1, Pa2, Pa4
}

M_matrix = [
    [0, 1, 0, 0, 1, 0],  # Pa1 -> Pa2, Pa5
    [0, 0, 1, 0, 0, 1],  # Pa2 -> Pa3, Pa6
    [0, 0, 0, 0, 0, 0],  # Pa3
    [0, 0, 0, 0, 0, 0],  # Pa4
    [0, 0, 0, 0, 0, 0],  # Pa5
    [1, 0, 0, 0, 0, 0]   # Pa6 -> Pa1
]

Pa_names = ['Pa1', 'Pa2', 'Pa3', 'Pa4', 'Pa5', 'Pa6']

Pa_labels = {
    'Pa1': 'Wind-Cold Pathogen (風寒病邪)',
    'Pa2': 'Stasis Blood Pathogen (瘀血病邪)',
    'Pa3': 'Lung Qi Deficiency (肺氣虛)',
    'Pa4': 'Lung Heat Pathogen (肺熱病邪)',
    'Pa5': 'Stomach Heat Pathogen (胃熱病邪)',
    'Pa6': 'Heart Qi Deficiency (心氣虛)'
}

#----------------------
# 2. 實作三個核心函式
#----------------------

# (1) 生成所有可能的病因網絡
def generate_possible_networks(symptoms, C_matrix, M_matrix):
    pa_sets = []
    for symptom in symptoms:
        # 找出 symptom 對應的 Pa
        pa_set = [Pa_names[i] for i, val in enumerate(C_matrix[symptom]) if val == 1]
        pa_sets.append(pa_set)
    
    # 將所有 pa_sets 做笛卡兒積，得到所有組合
    vertex_combinations = list(itertools.product(*pa_sets))
    
    PNS = []
    for combo in vertex_combinations:
        V_i = list(set(combo))  # 頂點集合 (去重複)
        
        # 根據 M_matrix 生成邊
        E_i = []
        for i, pa_i in enumerate(V_i):
            for j, pa_j in enumerate(V_i):
                if i != j:
                    idx_i = Pa_names.index(pa_i)
                    idx_j = Pa_names.index(pa_j)
                    if M_matrix[idx_i][idx_j] == 1:
                        E_i.append((pa_i, pa_j))
        
        # 組裝成一個網路資料結構
        G_i = {'vertices': V_i, 'edges': E_i}
        PNS.append(G_i)
    
    return PNS

# (2) 篩選有效病因網絡 (h-PNs)
def screen_valid_networks(PNS):
    h_PNS = []
    for G_i in PNS:
        G = nx.DiGraph()
        G.add_nodes_from(G_i['vertices'])
        G.add_edges_from(G_i['edges'])
        
        in_degrees = dict(G.in_degree())
        # 計算入度為 0 的節點個數
        N_G = sum(1 for v in in_degrees.values() if v == 0)
        
        if N_G >= 2:
            # Case 1: nh-PN
            continue
        
        elif N_G == 1:
            # Case 2: 入度為 0 的節點
            v_m = [v for v, deg in in_degrees.items() if deg == 0][0]
            reachable = set(nx.descendants(G, v_m)).union({v_m})
            if len(reachable) == len(G.nodes):
                G_i['root_vertex'] = v_m
                h_PNS.append(G_i)
        
        elif N_G == 0:
            # Case 3: 檢查是否有環 + 環中節點能到達所有節點
            try:
                cycle = nx.find_cycle(G)
                v_m = cycle[0][0]
                reachable = set(nx.descendants(G, v_m)).union({v_m})
                if len(reachable) == len(G.nodes):
                    G_i['root_vertex'] = v_m
                    h_PNS.append(G_i)
            except nx.NetworkXNoCycle:
                continue
    return h_PNS

# (3) 選擇頂點最少的網絡作為最終網絡
def select_final_network(h_PNS):
    if not h_PNS:
        return None
    return min(h_PNS, key=lambda x: len(x['vertices']))

#----------------------
# 3. 繪製病機演變圖
#----------------------
def draw_pathogenesis_network(network, symptoms):
    """
    network: {'vertices': [...], 'edges': [...], 'root_vertex': ...}
    symptoms: 輸入症狀清單，用於畫出紅色節點與邊
    """
    if not network:
        print("無有效病因網絡可繪製")
        return
    
    G = nx.DiGraph()
    
    # 症狀節點（紅色）
    symptom_nodes = [(s, {'color': 'red'}) for s in symptoms]
    # 病因節點（藍色，使用 Pa_labels 做可讀名稱）
    pa_nodes = [(Pa_labels[pa], {'color': 'blue'}) for pa in network['vertices']]
    
    G.add_nodes_from(symptom_nodes)
    G.add_nodes_from(pa_nodes)
    
    # 症狀 -> 病因原子
    S_to_Pa_edges = []
    for symptom in symptoms:
        for i, val in enumerate(C_matrix[symptom]):
            if val == 1:
                S_to_Pa_edges.append((symptom, Pa_labels[Pa_names[i]]))
    G.add_edges_from(S_to_Pa_edges)
    
    # 病因原子間的因果關係 (Pa -> Pa)
    G.add_edges_from([
        (Pa_labels[pa1], Pa_labels[pa2]) for (pa1, pa2) in network['edges']
    ])
    
    # 取出每個節點的 color 屬性
    node_colors = [data['color'] for node, data in G.nodes(data=True)]
    
    pos = nx.spring_layout(G, seed=42)  # 固定布局
    nx.draw(G, pos,
            with_labels=True,
            node_color=node_colors,
            node_size=2000,
            font_size=10,
            font_weight='bold',
            edge_color='gray',
            arrows=True,
            arrowstyle='->',
            arrowsize=20)
    
    # 設置中文字體避免警告
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
    plt.rcParams['axes.unicode_minus'] = False
    
    plt.title("Pathogenesis Evolution Graph (病機演變圖)")
    plt.show()

#----------------------
# 4. 主程式：執行並繪圖
#----------------------
if __name__ == "__main__":
    # 輸入症狀
    input_symptoms = ['S1', 'S2', 'S3']
    
    # 1) 生成所有可能網絡
    PNS = generate_possible_networks(input_symptoms, C_matrix, M_matrix)
    # 2) 篩選有效網絡
    h_PNS = screen_valid_networks(PNS)
    # 3) 選擇最終網絡
    final_network = select_final_network(h_PNS)
    
    print(f"Possible Networks: {len(PNS)}")
    print(f"Valid Networks (h-PNs): {len(h_PNS)}")
    if final_network:
        print(f"Final Network Vertices: {final_network['vertices']}")
        print(f"Final Network Edges: {final_network['edges']}")
        print(f"Root Vertex: {final_network.get('root_vertex', None)}")
    else:
        print("No valid network found.")
    
    # 4) 繪製圖
    draw_pathogenesis_network(final_network, input_symptoms)
