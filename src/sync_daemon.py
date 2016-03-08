#!/usr/bin/env python3

#Source:
# http://www.jejik.com/articles/2007/02/a_simple_unix_linux_daemon_in_python/
# Download in March 8th, 2016

#With custom modifications by
#Ataias Pereira Reis

import sys, time
from daemon import Daemon

import subprocess

class MyDaemon(Daemon):
    def run(self):
        while(True):
            time.sleep(30)
            subprocess.call(["/Users/ataias/Documents/ferrofluidos/src/dropbox_sync.py", "data", "/Users/ataias/Documents/ferrofluidos/src/data", "--token", "JA-2YSQygh4AAAAAAACKR7GGYzKj26F4TKBOtQOhE4XH_pz3DE7SEmr3R5mlElgr", "-y"])

if __name__ == "__main__":
    daemon = MyDaemon('/tmp/daemon-sync-daemon.pid')
    if len(sys.argv) == 2:
        if 'start' == sys.argv[1]:
            daemon.start()
        elif 'stop' == sys.argv[1]:
            daemon.stop()
        elif 'restart' == sys.argv[1]:
            daemon.restart()
        else:
            print("Unknown command")
            sys.exit(2)
        sys.exit(0)
    else:
        print("usage: %s start|stop|restart" % sys.argv[0])
        sys.exit(2)
