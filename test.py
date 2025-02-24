import itertools
import networkx as nx
import matplotlib.pyplot as plt

# 更新後的C矩陣
C_matrix = {
    'S1': [1, 0, 0, 0, 1, 0],  # S1 -> Pa1, Pa5 (e.g., 發燒伴寒戰 -> 風寒病邪, 胃熱病邪)
    'S2': [1, 0, 0, 0, 0, 1],  # S2 -> Pa1, Pa6 (e.g., 頭痛 -> 風寒病邪, 心氣虛)
    'S3': [1, 1, 0, 0, 0, 0],  # S3 -> Pa1, Pa2 (e.g., 舌淡 -> 風寒病邪, 瘀血病邪)
    'S4': [1, 0, 1, 1, 0, 0],  # S4 -> Pa1, Pa3, Pa4 (e.g., 浮脈 -> 風寒病邪, 肺氣虛, 肺熱病邪)
    'S5': [1, 0, 1, 1, 0, 0],  # S5 -> Pa1, Pa3, Pa4 (e.g., 薄白舌苔 -> 風寒病邪, 肺氣虛, 肺熱病邪)
    'S6': [1, 1, 0, 1, 0, 0]   # S6 -> Pa1, Pa2, Pa4 (e.g., 多汗 -> 風寒病邪, 瘀血病邪, 肺熱病邪)
}

# M矩陣保持不變
M_matrix = [
    [0, 1, 0, 0, 1, 0],  # Pa1 -> Pa2, Pa5 (風寒病邪 -> 瘀血病邪, 胃熱病邪)
    [0, 0, 1, 0, 0, 1],  # Pa2 -> Pa3, Pa6 (瘀血病邪 -> 肺氣虛, 心氣虛)
    [0, 0, 0, 0, 0, 0],  # Pa3 (肺氣虛)
    [0, 0, 0, 0, 0, 0],  # Pa4 (肺熱病邪)
    [0, 0, 0, 0, 0, 0],  # Pa5 (胃熱病邪)
    [1, 0, 0, 0, 0, 0]   # Pa6 -> Pa1 (心氣虛 -> 風寒病邪)
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

# 輸入症狀
input_symptoms = ['S1', 'S2', 'S3']

# 第一階段：生成所有可能的病因網絡
def generate_possible_networks(symptoms, C_matrix, M_matrix):
    pa_sets = []
    for symptom in symptoms:
        pa_set = [Pa_names[i] for i, val in enumerate(C_matrix[symptom]) if val == 1]
        pa_sets.append(pa_set)
    
    vertex_combinations = list(itertools.product(*pa_sets))
    PNS = []
    for combo in vertex_combinations:
        V_i = list(set(combo))
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

# 第二階段：Algorithm 3 - 篩選有效病因網絡 (h-PNs)
def screen_valid_networks(PNS):
    h_PNS = []
    
    for G_i in PNS:
        G = nx.DiGraph()
        G.add_nodes_from(G_i['vertices'])
        G.add_edges_from(G_i['edges'])
        
        in_degrees = dict(G.in_degree())
        N_G = sum(1 for v in in_degrees.values() if v == 0)
        
        if N_G >= 2:  # Case 1: nh-PN
            continue
        
        elif N_G == 1:  # Case 2: 檢查入度為0的頂點
            v_m = [v for v, deg in in_degrees.items() if deg == 0][0]
            reachable = set(nx.descendants(G, v_m)).union({v_m})
            if len(reachable) == len(G.nodes):
                G_i['root_vertex'] = v_m
                h_PNS.append(G_i)
        
        elif N_G == 0:  # Case 3: 檢查環中的頂點
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

# 第三階段：選擇最少Pa的網絡
def select_final_network(h_PNS):
    if not h_PNS:
        return None
    return min(h_PNS, key=lambda x: len(x['vertices']))

# 繪製病機演變圖（修正後）
def draw_pathogenesis_network(network):
    if not network:
        print("無有效病因網絡可繪製")
        return
    
    # 創建有向圖
    G = nx.DiGraph()
    
    # 添加節點：明確指定每個節點的color屬性
    symptom_nodes = [(s, {'color': 'red'}) for s in input_symptoms]  # 症狀用紅色
    pa_nodes = [(Pa_labels[pa], {'color': 'blue'}) for pa in network['vertices']]  # 病因原子用藍色
    G.add_nodes_from(symptom_nodes)
    G.add_nodes_from(pa_nodes)
    
    # 添加邊：症狀到病因原子的映射（根據C_matrix）
    S_to_Pa_edges = []
    for symptom in input_symptoms:
        for i, val in enumerate(C_matrix[symptom]):
            if val == 1:
                S_to_Pa_edges.append((symptom, Pa_labels[Pa_names[i]]))
    G.add_edges_from(S_to_Pa_edges)
    
    # 添加邊：病因原子間的因果關係
    G.add_edges_from([(Pa_labels[pa1], Pa_labels[pa2]) for pa1, pa2 in network['edges']])
    
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

# 中醫辨證系統（用於保留指定輸出）
def tcm_diagnosis_system(symptoms):
    # 假設的辨證邏輯
    if 'S1' in symptoms and 'S2' in symptoms:
        syndrome = "Wind-Cold Syndrome (風寒證)"
    elif 'S3' in symptoms or 'S4' in symptoms:
        syndrome = "Blood Stasis or Lung-Related Syndrome (瘀血或肺系證)"
    else:
        syndrome = "Unknown Syndrome (未知證候)"
    
    # 第一階段
    PNS = generate_possible_networks(symptoms, C_matrix, M_matrix)
    print(f"生成 {len(PNS)} 個可能的病因網絡")
    
    # 第二階段
    h_PNS = screen_valid_networks(PNS)
    print(f"篩選出 {len(h_PNS)} 個有效病因網絡")
    
    # 第三階段
    final_network = select_final_network(h_PNS)
    
    # 輸出
    print(f"輸入症狀: {symptoms}")
    print(f"辨證結果: {syndrome}")
    if final_network:
        print(f"最終病因網絡: Vertices: {final_network['vertices']}, Edges: {final_network['edges']}, Root Vertex: {final_network['root_vertex']}")
    else:
        print("無有效病因網絡")
    
    # 繪製病機圖
    draw_pathogenesis_network(final_network)

# 測試，保留指定輸出
if __name__ == "__main__":
    input_symptoms = ['S1', 'S2', 'S3']
    tcm_diagnosis_system(input_symptoms)