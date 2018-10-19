#!/usr/bin/env python

# File Send/Receive

from __future__ import print_function
import os
from sys import argv

if __name__ != "__main__":
    raise RuntimeError('Must run this file as the main program!')
if "TMP_HOST" not in os.environ:
    raise RuntimeError("must set TMP_HOST environ variable!")

def canwrite(xdir):
    uid, gid = os.getuid(), os.getgid()
    fstat = os.stat(xdir)
    fuid, fgid = fstat.st_uid, fstat.st_gid
    fmode = fstat.st_mode
    if fuid == uid:
        mask = 0o200
    elif fgid == gid:
        mask = 0o020
    else:
        mask = 0o002
    return (fmode & mask) != 0

def send_recv(args, task):
    num = "0"
    pnum = "./"
    ctmp = False
    dtmp = "%s/TMP" % os.environ["HOME"]
    if len(args) != 0 and args[0].isdigit():
        num = args[0]
        args = args[1:]
    assert len(args) == 0 or task == "send"
    if not canwrite("."):
        if not os.path.exists(dtmp):
            os.mkdir(dtmp)
            ctmp = True
        pnum = dtmp + "/"
    if task == "send":
        cmd = "zip %s.zip %s" % (pnum + num, " ".join(args))
        print (os.popen(cmd).read().strip())
        cmd = "scp %s.zip %s" % (pnum + num, os.environ["TMP_HOST"])
        for l in os.popen(cmd):
            print (l.strip())
    else:
        cmd = "scp %s/%s.zip %s" % (os.environ["TMP_HOST"], num, pnum)
        for l in os.popen(cmd):
            print (l.strip())
        cmd = "unzip %s.zip" % (pnum + num)
        print (os.popen(cmd).read().strip())
    cmd = "rm %s.zip" % (pnum + num)
    print (os.popen(cmd).read().strip())
    if ctmp:
        print (os.popen("rm -r " + dtmp).read().strip())

if argv[0].endswith("zs"):
    send_recv(argv[1:], "send")
elif argv[0].endswith("zr"):
    send_recv(argv[1:], "recv")
else:
    raise RuntimeError("This script cannot be run directly!")
