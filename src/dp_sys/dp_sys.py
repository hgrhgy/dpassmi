import numpy as np


class DoublePendulumSystem:
    """
    双摆锤系统定义
    """
    def __init__(self, l1, l2, m1, m2):
        """
        构造函数
        :param l1: 锤摆1长度
        :param l2: 锤摆2长度
        :param m1: 摆锤1质量
        :param m2: 摆锤2质量
        """
        self.l1 = l1
        self.l2 = l2
        self.m1 = m1
        self.m2 = m2
