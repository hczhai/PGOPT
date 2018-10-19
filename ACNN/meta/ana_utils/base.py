
import os

# construct file structure in a zip file
class FileItem(object):
    def __init__(self, path):
        if path[-1] == '/':
            self.isfile = False
        else: self.isfile = True
        self.path = path
        if self.isfile:
            self.name = os.path.basename(self.path)
            self.dirname = os.path.dirname(self.path)
        else:
            self.name = os.path.basename(self.path[:-1])
            self.dirname = os.path.dirname(self.path[:-1])
        self.parents = [g for g in self.dirname.split('/') if len(g) != 0]
        self.subs = {}

    def path_names(self):
        return [sub.name for sub in self.subs.values() if not sub.isfile]

    def file_names(self):
        return [sub.name for sub in self.subs.values() if sub.isfile]

    def __repr__(self):
        if len(self.subs) != 0:
            return "{}({}): {}".format("File" if self.isfile else "Path", self.name,
                                       self.subs.values())
        else: return "{}({})".format("File" if self.isfile else "Path", self.name)

    def __getitem__(self, i):
        return self.subs[i]

    def __setitem__(self, i, v):
        self.subs[i] = v

# time representation
def time_span_short(time, hour=True):
    time = int(time)
    xsec = time % 60
    xmin = (time / 60) % 60
    xhours = time / 60 / 60
    if not hour:
        xmin += xhours * 60
    rstr = []
    rstr += [xsec]
    if xmin != 0 or xhours != 0:
        rstr += [xmin]
    if xhours != 0 and hour:
        rstr += [xhours]
    rstr = ["%02d" % x for x in rstr[:-1]] + ["%d" % rstr[-1]]
    return ':'.join(rstr[::-1])

def time_span_str(time, no_sec=False):
    time = int(time)
    xsec = time % 60
    xmin = (time / 60) % 60
    xhours = (time / 60 / 60) % 24
    xdays = time / 60 / 60 / 24
    rstr = []
    if not no_sec:
        rstr.append("%d seconds" % xsec)
    if xmin != 0:
        rstr.append("%d minutes " % xmin)
    if xhours != 0:
        rstr.append("%d hours " % xhours)
    if xdays != 0:
        rstr.append("%d days " % xdays)
    return ' '.join(rstr[::-1])

# time intervals and their concatenation
class TimeItem(object):
    def __init__(self):
        self.intervals = []

    def sort(self):
        for i in range(len(self.intervals) - 1):
            if i >= len(self.intervals) - 1:
                break
            if self.intervals[i+1][0] > self.intervals[i][1]:
                continue
            elif self.intervals[i+1][0] <= self.intervals[i][1]:
                if self.intervals[i+1][1] <= self.intervals[i][1]:
                    self.intervals = self.intervals[:i+1] + self.intervals[i+2:]
                else:
                    self.intervals[i][1] = self.intervals[i+1][1]
                    self.intervals = self.intervals[:i+1] + self.intervals[i+2:]

    def add(self, start, end):
        ind = 0
        while ind != len(self.intervals):
            if self.intervals[ind][0] > start:
                break
            ind += 1
        self.intervals.insert(ind, [start, end])
        self.sort()

    def has(self, t):
        for s, e in self.intervals:
            if t >= s and t < e:
                return True
        return False

    def adds(self, xs):
        for start, end in xs:
            self.add(start, end)

    def total_time(self):
        return self.intervals[-1][1] - self.intervals[0][0]

    def start(self):
        return self.intervals[0][0]

    def end(self):
        return self.intervals[-1][1]

    def running_time(self):
        return sum([x[1] - x[0] for x in self.intervals])

    def idle_time(self):
        return self.total_time() - self.running_time()

# time and step trajectory of each structure
class StepAndTime(object):
    def __init__(self, tid, mstep, mac):
        self.traj = []
        self.tid = tid
        self.mstep = mstep
        self.mac = mac

    def add(self, step, timex):
        self.traj.append([step, timex])

    def total_step(self):
        return sum([x[0] for x in self.traj])

    def total_time(self):
        return sum([x[1] for x in self.traj])

class DictItem(dict):
    def __getattr__(self, attrname):
        return self[attrname]

    def __setattr__(self, attrname, value):
        self[attrname] = value

    def __delattr__(self, attrname):
        del self[attrname]

class RawTimeRecord(object):
    Tid = -1
    Rstep = 0
    Mstep = -1
    Mac = False # accpted mc step
    Rtime = 0.0 # cost relax time
    Stime = 0.0 # start time
    def __init__(self, tid, rstep, mstep, mac, rtime, stime):
        self.Tid, self.Rstep, self.Mstep, self.Rtime = tid, rstep, mstep, rtime
        self.Stime, self.Mac = stime, mac
