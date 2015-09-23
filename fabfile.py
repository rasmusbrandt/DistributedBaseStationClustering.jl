from fabric.api import *

# Forward SSH keys for github
env.forward_agent = True
env.user = 'rabr5411'
env.code_dir = '/home/rabr5411/sims/DistributedBaseStationClustering.jl'
env.git_url = 'git@gitr.sys.kth.se:rabr5411/DistributedBaseStationClustering.jl.git'

def gastown():
    env.hosts = 'gastown.156106636.members.btmm.icloud.com'
    env.user = 'rasmus'
    env.code_dir = '/Users/rasmus/Desktop/sims/DistributedBaseStationClustering.jl'

def KTHdebsim():
    env.hosts = '130.237.50.52'
    env.user = 'rabrdeb'
    env.port = 1024
    env.code_dir = '/home/rabrdeb/sims/DistributedBaseStationClustering.jl'

def sim401():
    env.hosts = 'sim401.ee.kth.se'

def sim402():
    env.hosts = 'sim402.ee.kth.se'

def sim403():
    env.hosts = 'sim403.ee.kth.se'

def sim404():
    env.hosts = 'sim404.ee.kth.se'

def sim405():
    env.hosts = 'sim405.ee.kth.se'

def all():
    env.hosts = [ 'sim40{}.ee.kth.se'.format(id) for id in range(1,6) ]

def deploy(ref='master'):
    with settings(warn_only = True):
        if run('test -d {}'.format(env.code_dir)).failed:
            run('mkdir -p {}'.format(env.code_dir))
            run('git clone {0} {1}'.format(env.git_url, env.code_dir))
    with cd(env.code_dir):
        run('git fetch')
        run('git checkout {}'.format(ref))
        run('git pull')
