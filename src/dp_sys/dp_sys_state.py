from src.dp_sys.utils import deg2rad


class DoublePendulumState:

    def __init__(self, t, th1, th2, th1_rate, th2_rate):
        """
        双摆锤状态定义类
        :param t: 状态时刻
        :param th1: 锤摆1弧度
        :param th2: 锤摆2弧度
        :param th1_rate: 锤摆1角速度
        :param th2_rate: 锤摆2角速度
        """

        self.t = t
        self.th1 = th1
        self.th2 = th2
        self.th1_rate = th1_rate
        self.th2_rate = th2_rate


