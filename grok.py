import itertools
import networkx as nx
import matplotlib.pyplot as plt

# 定義 C矩陣(Sy2Pa)和 M矩陣(Pa2Pa)
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

# 第一階段：生成所有可能的病因網絡
def generate_possible_networks(symptoms, C_matrix, M_matrix):
    pa_sets = []
    for symptom in symptoms:
        if symptom in C_matrix:
            pa_set = [Pa_names[i] for i, val in enumerate(C_matrix[symptom]) if val == 1]
            pa_sets.append(pa_set)
        else:
            pa_sets.append([])  # 如果症狀不在 C_matrix 中，跳過
    if not pa_sets:  # 如果沒有有效的症狀，直接返回空列表
        return []
    
    vertex_combinations = list(itertools.product(*pa_sets))
    PNS = []
    for combo in vertex_combinations:
        V_i = list(set(combo))  # 去重形成頂點集合
        E_i = []
        for i, pa_i in enumerate(V_i):
            for j, pa_j in enumerate(V_i):
                if i != j:
                    idx_i = Pa_names.index(pa_i)
                    idx_j = Pa_names.index(pa_j)
                    if M_matrix[idx_i][idx_j] == 1:
                        E_i.append((pa_i, pa_j))
        G_i = {'vertices': V_i, 'edges': E_i}
        PNS.append(G_i)
    return PNS

# 第二階段：Algorithm 3 - 篩選有效病因網絡 (h-PNs，修正後版本)
def screen_valid_networks(PNS):
    h_PNS = []
    
    for G_i in PNS:
        if not G_i['vertices']:  # 如果頂點集合為空，跳過
            continue
        
        G = nx.DiGraph()
        G.add_nodes_from(G_i['vertices'])
        G.add_edges_from(G_i['edges'])
        
        # 計算每個頂點的入度
        in_degrees = dict(G.in_degree())
        N_G = sum(1 for v in in_degrees.values() if v == 0)  # 入度為 0 的頂點數
        
        # Case 1: N_G >= 2, 根據 Theorem 1 標記為 nh-PN
        if N_G >= 2:
            continue
        
        # Case 2: N_G == 1, 檢查入度為 0 的頂點是否為 Rv
        elif N_G == 1:
            v_m = [v for v, deg in in_degrees.items() if deg == 0][0]
            reachable = set(nx.descendants(G, v_m)).union({v_m})
            if len(reachable) == len(G.nodes):  # v_m 是 Rv
                G_i['root_vertex'] = v_m
                h_PNS.append(G_i)
        
        # Case 3: N_G == 0, 檢查環中的頂點是否為 Rv
        elif N_G == 0:
            try:
                cycle = nx.find_cycle(G)
                v_m = cycle[0][0]  # 從環中選一個頂點
                reachable = set(nx.descendants(G, v_m)).union({v_m})
                if len(reachable) == len(G.nodes):  # v_m 是 Rv
                    G_i['root_vertex'] = v_m
                    h_PNS.append(G_i)
            except nx.NetworkXNoCycle:
                continue
    
    return h_PNS

# 第三階段：選擇邊數最多的有效網絡（修正後）
def select_final_network(h_PNS):
    if not h_PNS:
        return None
    return max(h_PNS, key=lambda x: len(x['edges']))  # 選擇邊最多的網絡

# 繪製病機演變圖（修正後，確保所有節點有 color 屬性）
def draw_pathogenesis_network(network, input_symptoms):
    if not network:
        print("無有效病因網絡可繪製")
        return
    
    # 創建有向圖
    G = nx.DiGraph()
    
    # 確保添加所有節點並設置 color 屬性
    symptom_nodes = [(s, {'color': 'red'}) for s in input_symptoms if s in C_matrix]  # 只有有效的症狀
    pa_nodes = [(Pa_labels[pa], {'color': 'blue'}) for pa in network['vertices'] if pa in Pa_names]  # 只有有效的 Pa
    G.add_nodes_from(symptom_nodes)
    G.add_nodes_from(pa_nodes)
    
    # 添加邊：症狀到病因原子的映射（根據 C_matrix）
    S_to_Pa_edges = []
    for symptom in input_symptoms:
        if symptom in C_matrix:
            for i, val in enumerate(C_matrix[symptom]):
                if val == 1:
                    pa = Pa_names[i]
                    if pa in network['vertices']:  # 僅添加映射到最終網絡中的 Pa
                        S_to_Pa_edges.append((symptom, Pa_labels[pa]))
    G.add_edges_from(S_to_Pa_edges)
    
    # 添加邊：病因原子間的因果關係
    G.add_edges_from([(Pa_labels[pa1], Pa_labels[pa2]) for pa1, pa2 in network['edges'] 
                      if pa1 in network['vertices'] and pa2 in network['vertices']])
    
    # 檢查並設置節點的 color 屬性，確保所有節點都有 color
    for node, data in G.nodes(data=True):
        if 'color' not in data:
            # 根據節點名稱判斷是症狀還是 Pa
            if node in input_symptoms:
                data['color'] = 'red'
            elif any(node.startswith(label) for label in Pa_labels.values()):
                data['color'] = 'blue'
    
    # 設置節點樣式
    node_colors = [data['color'] for node, data in G.nodes(data=True)]
    
    # 繪製圖形
    pos = nx.spring_layout(G, seed=42)  # 固定布局以確保一致
    nx.draw(G, pos, with_labels=True, node_color=node_colors, 
            node_size=2000, font_size=12, font_weight='bold', 
            edge_color='gray', arrows=True, arrowstyle='->', arrowsize=20)
    
    # 設置中文字體以避免警告
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']  # 支援中文
    plt.rcParams['axes.unicode_minus'] = False
    
    plt.title("Pathogenesis Evolution Graph (病機演變圖)")
    plt.show()

# 主程式：允許用戶輸入症狀並輸出有效病機網絡
def main():
    print("可輸入的症狀：S1, S2, S3, S4, S5, S6（用逗號分隔，例如 'S1,S3,S4'）")
    symptoms_input = input("請輸入症狀：").strip().split(',')
    input_symptoms = [s.strip() for s in symptoms_input if s.strip() in C_matrix]
    
    if not input_symptoms:
        print("無效的症狀輸入，請確保輸入 S1 到 S6 之間的症狀")
        return
    
    # 第一階段
    PNS = generate_possible_networks(input_symptoms, C_matrix, M_matrix)
    print(f"生成 {len(PNS)} 個可能的病因網絡")
    
    # 第二階段
    h_PNS = screen_valid_networks(PNS)
    print(f"篩選出 {len(h_PNS)} 個有效病因網絡")
    
    # 第三階段
    final_network = select_final_network(h_PNS)
    
    # 輸出結果
    print(f"輸入症狀: {input_symptoms}")
    if final_network:
        print(f"最終病因網絡: Vertices: {final_network['vertices']}, Edges: {final_network['edges']}, Root Vertex: {final_network['root_vertex']}")
    else:
        print("無有效病因網絡")
    
    # 繪製病機圖
    draw_pathogenesis_network(final_network, input_symptoms)

if __name__ == "__main__":
    main()