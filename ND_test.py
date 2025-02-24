import itertools

#C, M矩陣
C_matrix = {
    'S1': [1, 0, 0, 0, 1, 0],
    'S2': [1, 0, 0, 0, 0, 1],
    'S3': [1, 1, 0, 0, 0, 0],
    'S4': [1, 0, 1, 1, 0, 0],
    'S5': [1, 0, 1, 1, 0, 0],
    'S6': [1, 1, 0, 1, 0, 0]
}

M_matrix = [
    [0, 1, 0, 0, 1, 0],
    [0, 0, 1, 0, 0, 1],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 0] 
]

Pa_names = ['Pa1', 'Pa2', 'Pa3', 'Pa4', 'Pa5', 'Pa6']
input_symptoms = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']

def generate_possible_networks(symptoms, C_matrix, M_matrix):
    #step 1: 從症狀中找到映射病機原子
    pa_sets = []
    for symptom in symptoms:
        pa_set = [Pa_names[i] for i, val in enumerate(C_matrix[symptom]) if val==1]
        pa_sets.append(pa_set)

    #笛卡爾積生成所有可能的Pa組合
    vertex_combinations = list(itertools.product(*pa_sets))

    PNS = []
    for combo in vertex_combinations:
        V_i = list(set(combo)) #頂點集合

        #step 2 生成邊集合
        E_i = []
        for i, pa_i in enumerate(V_i):
            for j, pa_j in enumerate(V_i):
                if i != j:
                    idx_i = Pa_names.index(pa_i)
                    idx_j = Pa_names.index(pa_j)
                    if M_matrix[idx_i][idx_j] == 1:
                        E_i.append((pa_i, pa_j))
        
        #step 3 生成ND
        G_i = {'vertices': V_i, 'edges': E_i}
        PNS.append(G_i)

    return PNS

PNS = generate_possible_networks(input_symptoms, C_matrix, M_matrix)
for i, network in enumerate(PNS):
    print(f"G{i+1}:")
    print(f"  Vertices: {network['vertices']}")
    print(f"  Edges: {network['edges']}")