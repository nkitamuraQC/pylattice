import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from sklearn.linear_model import LinearRegression

class MD_Properties:
    def __init__(self, md_analyzers_list):
        self.md_analyzers = md_analyzers_list
        self.target_data = []
    
    
    def gather_data(self, target="stresses_au", index=None):
        if index is not np.array:
            raise NotImplementedError
        for md_analyser in self.md_analyzers:
            val = md_analyser[target][-1][index]
            self.target_data.append(val)
        self.target_data = np.array(self.target_data)
        return
    
    def plot(self):
        steps = [i for i in range(self.target_data)]
        plt.plot(steps, self.target_data)
        plt.savefig("plot.png")
        return


    def get_young(self, dlat, targets_list=None, target="stresses_au", index=None):
        self.gather_data(target=target, index=index)
        if len(self.target_data) == 0:
            raise ValueError
        elif targets_list is list:
            targets_array = np.array(targets_list)
            targets = self.target_data[targets_array]
            ds = dlat[index] * targets_array
            lr = LinearRegression()
            lr.fit(ds, targets)
            print('coefficient = ', lr.coef_[0]) # 説明変数の係数を出力
            print('intercept = ', lr.intercept_) 
            young = lr.coef_[0]
            return young

        return
    


    


