from seekpath import get_path

class SeekPath:
    def __init__(self, qe_ctrl):
        self.qe_ctrl = qe_ctrl
        self.coord = None
        self.path = None

    def exec(self):
        cell = self.qe_ctrl.lattice
        pos = self.qe_ctrl.xyz
        n = [i for i in range(len(pos))]
        struct = (cell, pos, n)
        res = get_path(struct)
        self.coord = res["point_coords"]
        self.path = res['path']
        
        path_sym = []
        for i, p in enumerate(self.path):
            path_sym.append(p[0])
            if i == len(self.path) - 1:
                path_sym.append(p[1])

        path_coord = []
        for p in path_sym:
            path_coord.append(self.coord[p])

        return path_sym, path_coord


