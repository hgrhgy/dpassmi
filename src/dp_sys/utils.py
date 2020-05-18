import math


def deg2rad(deg):
    """
    角度转弧度
    :param deg: 角度
    :return:
    """
    return deg * math.pi / 180.0


def rad2deg(rad):
    """
    弧度转角度
    :param rad: 弧度
    :return:
    """
    return rad * 180.0 / math.pi


def map2pi(rad):
    """
    弧度转化，把不在[-PI, PI]内的弧度转化为[-PI, PI]
    :param rad:
    :return:
    """
    if rad < -math.pi:
        return rad + 2 * math.pi
    elif rad > math.pi:
        return rad - 2 * math.pi
    return rad
