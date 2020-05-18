import matplotlib.pyplot as plt

def plot_states(states):
    plt.plot(range(0, len(states)), [s.th1 for s in states])
    plt.show()


def plot_states_and_observation(states, obs_states):
    plt.plot(range(0, len(states)), [s.th1 for s in states])
    plt.scatter([s.t for s in obs_states], [s.th1 for s in obs_states], color='r')
    plt.show()