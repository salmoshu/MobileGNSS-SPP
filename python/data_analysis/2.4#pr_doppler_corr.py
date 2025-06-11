import matplotlib.pyplot as plt

# 对应文档的2.4节（伪距和多普勒之间的相关性）

# 运行RTKLIB时，需要将钟差和钟漂信息打印到pos文件中（默认的pos格式是不包含钟差和钟漂信息的）
def get_clk_from_pos(pos_path):
    with open(pos_path, 'r') as f:
        lines = f.readlines()
        clk_bias = []
        clk_drift = []
        for line in lines:
            if line.startswith('%'):
                continue
            fields = line.split()
            clk_bias.append(float(fields[-2]))
            clk_drift.append(float(fields[-1]))

    return clk_bias, clk_drift

if __name__ == "__main__":
    path = r'./rover_clk.pos'
    bias, drift = get_clk_from_pos(path)

    idx = 0
    idx_arr = []
    bias_drift = [None]*len(bias)
    for b, d in zip(bias, drift):
        if idx == 0:
            bias_drift[idx] = None
        else:
            last_b = bias[idx-1]
            bias_drift[idx] = b - last_b

        idx_arr.append(idx)
        idx += 1
    
    plt.plot(idx_arr[1:], bias_drift[1:], label='bias_drift')
    plt.plot(idx_arr[1:], drift[1:], label='drift')
    plt.legend()
    plt.grid(True)
    plt.show()
