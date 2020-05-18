from src.dp_sys.dp_sys import DoublePendulumSystem
from src.dp_sys.dp_sys_obs import DoublePendulumObserver
from src.dp_sys.dp_sys_plot import plot_states, plot_states_and_observation
from src.dp_sys.dp_sys_sim import DoublePendulumSimulator
from src.dp_sys.dp_sys_state import DoublePendulumState
from src.dp_sys.utils import deg2rad
from src.dp_sys_config import *

if __name__ == '__main__':
    init_state = DoublePendulumState(0, deg2rad(angle1), deg2rad(angle2), deg2rad(roc_angle1), deg2rad(roc_angle2))
    system = DoublePendulumSystem(l1, l2, m1, m2)
    simulator = DoublePendulumSimulator(system, init_state, g, obs_window, length_of_integration, time_step)
    observer = DoublePendulumObserver(obs_freq, deg2rad(noise_pos), deg2rad(noise_vel))

    simulator.calc_all_states()
    observer.make_obs(simulator, 0, 2000)

    plot_states(simulator.states_in_duration)
    plot_states_and_observation(simulator.states_in_duration, observer.obs_states)

