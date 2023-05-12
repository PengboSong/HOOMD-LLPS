import datetime
from string import Template


class TimeDeltaTemplate(Template):
    delimiter = '%'


def strfdelta(tdelta: datetime.timedelta, fmt: str):
    d = {'D': tdelta.days}
    d['H'], tleft = divmod(tdelta.seconds, 3600)
    d['M'], d['S'] = divmod(tleft, 60)
    return TimeDeltaTemplate(fmt).substitute(**d)