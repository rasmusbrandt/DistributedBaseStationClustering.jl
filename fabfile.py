from fabric.api import *

env.forward_agent = True
env.git_url = 'git@gitr.sys.kth.se:rabr5411/IAClustering.jl.git'

def gastown():
    env.hosts = 'gastown.156106636.members.btmm.icloud.com'
    env.user = 'rasmus'
    env.code_dir = '/Users/rasmus/Desktop/sims/IAClustering.jl'

def deploy(ref='master'):
    with settings(warn_only = True):
        if run('test -d {}'.format(env.code_dir)).failed:
            run('mkdir -p {}'.format(env.code_dir))
            run('git clone {0} {1}'.format(env.git_url, env.code_dir))
    with cd(env.code_dir):
        run('git pull')
        run('git checkout {}'.format(ref))
