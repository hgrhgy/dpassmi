import numpy as np

from src.dp_sys.dp_sys_obs_state import DoublePendulumObservationState


class DoublePendulumObserver:

    def __init__(self, freq, noi_pos, noi_vel):
        """
        观测者
        :param freq: 观测频率
        :param noi_pos: 观测角度误差， 弧度
        :param noi_vel: 观测角速度误差， 弧度
        """
        self.freq = freq
        self.noi_pos = noi_pos
        self.noi_vel = noi_vel
        self.obs_states = []
        pass

    def make_obs(self, simulator, start_time, end_time):
        """
        观测
        :param simulator: 系统模拟器
        :param start_time: 观测起始时间
        :param end_time: 观测结束时间
        :param freq: 观测频率
        :return:
        """
        freq = self.freq
        noi_pos = self.noi_pos
        noi_vel = self.noi_vel

        sarr = []
        for i in range(start_time, end_time):
            if i % freq == 0:
                sarr.append(simulator.states_in_duration[i])

        for s in sarr:
            self.obs_states.append(DoublePendulumObservationState(s.t, self.observer_with_noise(s.th1, noi_pos),
                                                                  self.observer_with_noise(s.th2, noi_pos)))

    def observer_with_noise(self, v, n):
        """
        带误差观测
        :param v: 观测值
        :param n: 最大误差
        :return:
        """
        return v + (2 * np.random.rand() - 0.5) * n
