


class DoublePendulumObservationState:

    def __init__(self, t, th1, th2):
        """
        双摆锤观测状态定义类
        :param t: 状态时刻
        :param th1: 锤摆1弧度
        :param th2: 锤摆2弧度
        """

        self.t = t
        self.th1 = th1
        self.th2 = th2
