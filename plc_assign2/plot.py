import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# clock() version 1 
t = pd.read_csv("./log_cloud_4am.txt", header=1)
t1 = pd.read_csv("./log_cloud_ripple_4am.txt", header=1)
# clock() version 2 
t2 = pd.read_csv("./log_cloud.txt", header=1)
t3 = pd.read_csv("./log_cloud_ripple.txt", header=1)



# clock() version 1 

# plot 1: First, plot the execution time of 2, 4, 8, 16 and 32 rank runs as function of their number of ranks
# fig, axes = plt.subplots(3, 3, figsize=(10, 10))

# fig.delaxes(axes[-1][-1])
# fig.delaxes(axes[-1][-2])

# i = 0
# for test in range(1, 8):
#     axes[i//3][i%3].scatter([2, 4, 8, 16, 32], np.array(t[t.test==test][t.barrier==1].groupby(["rank"])["time"].mean()), color='darksalmon', label="with barrier")
#     axes[i//3][i%3].plot([2, 4, 8, 16, 32], np.array(t[t.test==test][t.barrier==1].groupby(["rank"])["time"].mean()), color='darksalmon')
    
#     axes[i//3][i%3].scatter([2, 4, 8, 16, 32], np.array(t[t.test==test][t.barrier==0].groupby(["rank"])["time"].mean()), color='firebrick', label="without barrier")
#     axes[i//3][i%3].plot([2, 4, 8, 16, 32], np.array(t[t.test==test][t.barrier==0].groupby(["rank"])["time"].mean()), color='firebrick')

#     axes[i//3][i%3].set_title("test=%s" % test)
#     axes[i//3][i%3].set_xticks([2, 4, 8, 16, 32])
#     axes[i//3][i%3].set_xticklabels(['2', '4', '8', '16', '32'])
#     axes[i//3][i%3].set_xlabel("rank")
#     axes[i//3][i%3].set_ylabel("time(seconds)")
#     axes[i//3][i%3].legend()

#     i += 1

# plt.tight_layout()

# fig.savefig("rank_time.png", dpi=300)


# plot 2: Next, plot of the speedup relative to the execution time of the serial MPI CLA adder of the 2, 4, 8, 16 and 32 rank runs

fig, axes = plt.subplots(3, 3, figsize=(10, 10))

fig.delaxes(axes[-1][-1])
fig.delaxes(axes[-1][-2])

i = 0
for test in range(1, 8):
    axes[i//3][i%3].scatter([2, 4, 8, 16, 32], np.array(t1[t.test==test]["time"].mean - 
                                                        t[t.test==test][t.barrier==1].groupby(["rank"])["time"].mean()
                                                            ), color='darksalmon', label="with barrier")
    # axes[i//3][i%3].plot([2, 4, 8, 16, 32], np.array(t[t.test==test][t.barrier==1].groupby(["rank"])["time"].mean()), color='darksalmon')

    # axes[i//3][i%3].scatter([2, 4, 8, 16, 32], np.array((t1[t.test==test][t.barrier==0].groupby(["rank"])["time"].mean() - 
    #                                                         t[t.test==test]["time"].mean()), color='darksalmon', label="without barrier")
    # axes[i//3][i%3].plot([2, 4, 8, 16, 32], np.array(t[t.test==test][t.barrier==0].groupby(["rank"])["time"].mean()), color='darksalmon')

    

    axes[i//3][i%3].set_title("test=%s" % test)
    axes[i//3][i%3].set_xticks([2, 4, 8, 16, 32])
    axes[i//3][i%3].set_xticklabels(['2', '4', '8', '16', '32'])
    axes[i//3][i%3].set_xlabel("rank")
    axes[i//3][i%3].set_ylabel("speedup( to serial MPI CLA adder)")
    axes[i//3][i%3].legend()

    i += 1

plt.tight_layout()

fig.savefig("rank_time.png", dpi=300)

# plot 3: plot the speedup of the relative to the execution time of the serial ripplecarry adder to the MPI CLA adder running in parallel on 2 thru 32 ranks.  