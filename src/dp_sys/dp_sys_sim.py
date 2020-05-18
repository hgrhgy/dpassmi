from src.dp_sys.dp_sys_math import runge_kutt_4
from src.dp_sys.dp_sys_state import DoublePendulumState


class DoublePendulumSimulator:

    def __init__(self, system, init_state, g, win, duration, dt):
        """
        双摆锤模拟器
        :param system: 双摆锤系统
        :param state: 系统初始状态
        :param g: 重力参数
        :param window: 窗口大小
        :param duration: 模拟总时长
        :param dt: 时间步长
        """
        self.system = system
        self.init_state = init_state
        self.g = g
        self.win = win
        self.duration = duration
        self.dt = dt

        # 根据窗口求模拟时长中总窗口数
        self.num_of_win = int(duration / win) + int(duration % win > 0)

        # 状态量按窗口大小划分
        self.states_by_win = []

        # 所有状态量
        self.states_in_duration = []

        # 当前状态,初始化时为初始状态
        self.state = init_state

    def calc_all_states(self):
        """
        计算所有状态
        :return:
        """
        s = self.init_state
        for i in range(0, self.duration):
            self.states_in_duration.append(s)
            tmp = runge_kutt_4(s.th1, s.th2, s.th1_rate, s.th2_rate, self.system.m1, self.system.m2, self.system.l1,
                                self.system.l2, self.g, self.dt)
            s = DoublePendulumState(i, tmp[0], tmp[1], tmp[2], tmp[3])
