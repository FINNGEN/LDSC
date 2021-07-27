import os,sys,subprocess

def make_sure_path_exists(path):
    import errno
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
def progressBar(value, endvalue, bar_length=20):
    '''
    Writes progress bar, given value (eg.current row) and endvalue(eg. total number of rows)
    '''

    percent = float(value) / endvalue
    arrow = '-' * int(round(percent * bar_length)-1) + '>'
    spaces = ' ' * (bar_length - len(arrow))

    sys.stdout.write("\rPercent: [{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
    sys.stdout.flush()

def tmp_bash(cmd,check = False):
    from tempfile import NamedTemporaryFile

    scriptFile = NamedTemporaryFile(delete=True)
    with open(scriptFile.name, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write(cmd + "\n")

    os.chmod(scriptFile.name,0o777)
    scriptFile.file.close()

    if check:
        subprocess.check_call(scriptFile.name)
    else:
        subprocess.call(scriptFile.name,stderr = subprocess.DEVNULL,stdout  = subprocess.DEVNULL)
