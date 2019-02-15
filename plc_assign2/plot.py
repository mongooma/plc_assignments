import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

t = pd.read_csv("./log_cloud_1.txt", header=1)

fig, axes = plt.subplots(3, 3, figsize=(10, 10))

fig.delaxes(axes[-1][-1])
fig.delaxes(axes[-1][-2])

i = 0
for test in range(1, 8):
    axes[i//3][i%3].scatter([2, 4, 8, 16, 32], np.array(t[t.test==test][t.barrier==1].groupby(["rank"])["time"].mean()), color='darksalmon', label="with barrier")
    axes[i//3][i%3].plot([2, 4, 8, 16, 32], np.array(t[t.test==test][t.barrier==1].groupby(["rank"])["time"].mean()), color='darksalmon')
    
    axes[i//3][i%3].scatter([2, 4, 8, 16, 32], np.array(t[t.test==test][t.barrier==0].groupby(["rank"])["time"].mean()), color='firebrick', label="without barrier")
    axes[i//3][i%3].plot([2, 4, 8, 16, 32], np.array(t[t.test==test][t.barrier==0].groupby(["rank"])["time"].mean()), color='darksalmon')

    axes[i//3][i%3].set_title("test=%s" % test)
    axes[i//3][i%3].set_xticks([2, 4, 8, 16, 32])
    axes[i//3][i%3].set_xticklabels(['2', '4', '8', '16', '32'])
    axes[i//3][i%3].set_xlabel("rank")
    axes[i//3][i%3].set_ylabel("time(seconds)")
    axes[i//3][i%3].legend()

    i += 1

plt.tight_layout()

fig.savefig("rank_time.png", dpi=300)