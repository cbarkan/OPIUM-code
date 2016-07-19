import numpy as np
from subprocess import Popen
import threading
filename = 'v'
#p = subprocess.Popen('./opium ' + filename + ' ' + filename +'.log all rpt', shell=True)

def run_command_with_timeout(cmd, timeout_sec):
    """Execute `cmd` in a subprocess and enforce timeout `timeout_sec` seconds.
 
    Return subprocess exit code on natural completion of the subprocess.
    Raise an exception if timeout expires before subprocess completes."""
    proc = Popen(cmd, shell=True)
    proc_thread = threading.Thread(target=proc.communicate)
    proc_thread.start()
    proc_thread.join(timeout_sec)
    if proc_thread.is_alive():
        # Process still running - kill it and raise timeout error
        try:
            proc.kill()
        except OSError, e:
            # The process finished between the `is_alive()` and `kill()`
            return proc.returncode
        # OK, the process was definitely killed
        #raise SubprocessTimeoutError('Process #%d killed after %f seconds' % (proc.pid, timeout_sec))
    # Process completed naturally - return exit code
    
cmd = './opium ' + filename + ' ' + filename +'.log all rpt'
run_command_with_timeout(cmd, 1)