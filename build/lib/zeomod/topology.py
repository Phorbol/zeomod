import numpy as np
from ase import Atoms
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
from matscipy.neighbours import neighbour_list

class GraphProcessor:
    """使用 matscipy 和 scipy.sparse 的高性能图处理器"""
    def __init__(self, atoms: Atoms):
        print(" [GraphProcessor] 初始化：使用 matscipy 构建T-T网络...")
        self.atoms = atoms
        self.t_indices = [a.index for a in self.atoms if a.symbol in ['Si', 'Al']]
        self.num_t_sites = len(self.t_indices)

        if self.num_t_sites == 0:
            self.distance_matrix = np.array([])
            return

        # 映射表
        self.idx_map = {original_idx: i for i, original_idx in enumerate(self.t_indices)}
        
        # 构建并计算
        graph_matrix = self._build_sparse_matrix()
        self.distance_matrix = shortest_path(csgraph=graph_matrix, directed=False, unweighted=True)
        self.distance_matrix[np.isinf(self.distance_matrix)] = -1
        self.distance_matrix = self.distance_matrix.astype(int)

    def _build_sparse_matrix(self) -> csr_matrix:
        """构建 T-O-T 连接矩阵"""
        # Matscipy 极速计算邻居
        i_idx, j_idx = neighbour_list('ij', self.atoms, cutoff=2.1) # 略微放大 cutoff 确保覆盖

        adj = [[] for _ in range(len(self.atoms))]
        for i, j in zip(i_idx, j_idx):
            adj[i].append(j)

        o_indices = [a.index for a in self.atoms if a.symbol == 'O']
        t_symbols = {'Si', 'Al'}

        row_ind, col_ind = [], []
        
        # 遍历所有氧原子寻找桥接
        for o_idx in o_indices:
            neighbors = adj[o_idx]
            t_neighbors = [idx for idx in neighbors if self.atoms[idx].symbol in t_symbols]
            
            if len(t_neighbors) == 2:
                t1 = self.idx_map[t_neighbors[0]]
                t2 = self.idx_map[t_neighbors[1]]
                row_ind.extend([t1, t2])
                col_ind.extend([t2, t1])

        data = np.ones(len(row_ind), dtype=np.int8)
        return csr_matrix((data, (row_ind, col_ind)), shape=(self.num_t_sites, self.num_t_sites))

    def get_distance(self, t_site1_idx: int, t_site2_idx: int) -> int:
        if t_site1_idx not in self.idx_map or t_site2_idx not in self.idx_map:
            return -1
        return self.distance_matrix[self.idx_map[t_site1_idx], self.idx_map[t_site2_idx]]