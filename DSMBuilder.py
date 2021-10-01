#
# SCons build script for building Fortran applications using the DSM API
#
# (c) 2021- Wim Vanderbauwhede <wim.vanderbauwhede@mail.be>
#


import os
import re
import subprocess
import sys
import os.path
from SCons.Variables import Variables
from SCons.Environment import Environment
opts=''

def getOpt(optname,desc,default):
    global opts
    opts.Add(optname,desc,default)
    optionl = list(filter (lambda x: x.key==optname,opts.options))
    if optionl:
        option=optionl[0]
        if optname in opts.args and opts.args[optname]!=option.default:
            print(opts.args)
            return opts.args[option.key]
        else: 
            return option.default
    else:
        print( "No such option: "+optname )
    
def initDSM(*envt):
    if envt==():
        env=Environment()
    else:
        env=envt[0]

    global opts, useF, useFC, useDyn, buildLib, installPrefix, mcModel, DSM_FORTRAN_API_DIR

    DSM_FORTRAN_API_DIR='.'
    if not 'DSM_FORTRAN_API_DIR' in os.environ:
        print('The environment variable DSM_FORTRAN_API_DIR must be defined')
        exit(1)
    else:
        DSM_FORTRAN_API_DIR=os.environ['DSM_FORTRAN_API_DIR']
    if 'DSM_FORTRAN_LIB_PREFIX' in os.environ:
        DSM_FORTRAN_LIB_PREFIX=os.environ['DSM_FORTRAN_LIB_PREFIX']
    else:
        DSM_FORTRAN_LIB_PREFIX=DSM_FORTRAN_API_DIR
 
    if subprocess.check_output(['uname',''], shell=True).strip()  == "Darwin":
        OSX=1
        OSFLAG='-DOSX'
    else:
        OSX=0
        OSFLAG='-D__LINUX__'

    opts=Variables()        
    CWD= os.environ['PWD']
    args=sys.argv[1:]
    for arg in args:
        if re.match("(\w+)=(\w+)",arg):
            (k,v)=arg.split('=')
            opts.args[k]=v
    env['KERNEL_OPTS']=[]   

    buildLib = getOpt('lib', 'Build and install the library','0')
    if 'buildLib' in env:
        if env['buildLib']==1:
            buildLib='1'

    installPrefix = DSM_FORTRAN_LIB_PREFIX
    if 'installPrefix' in env:
        installPrefix=env['installPrefix']
    
    if buildLib=='1':
        help = """
        Options:
        The environment variable DSM_FORTRAN_API_DIR must be defined as the path to the DSM Fortran API directory.
        If the environment variable DSM_FORTRAN_LIB_PREFIX is defined it will be used as the prefix for the library installation.
        It defaults to DSM_FORTRAN_API_DIR
        
        *dyn=0|1 [0] build a dynamic Library                            DSMBuilder.useDyn
        *mcm=s|m|l [s] mcmodel flag for gcc/gfortran                    DSMBuilder.mcModel
         O=[gfortran/gcc -O flag] [3]

        The options marked with * can be passed as arguments to the buildLibs call in your scons script
        The DSM Fortran API directory can be accessed via DSMBuilder.DSM_FORTRAN_API_DIR
        """ 
    else:
        help = """
        Options:
        The environment variable DSM_FORTRAN_API_DIR must be defined as the path to the DSM Fortran API directory.
        If the environment variable DSM_FORTRAN_LIB_PREFIX is defined it will be used as the prefix for the library installation.
        It defaults to DSM_FORTRAN_API_DIR

        *dyn=0|1 [0] build a dynamic Library             DSMBuilder.useDyn
        *mcm=s|m|l [s] mcmodel flag for gcc/gfortran
        O=[gcc -O flag] [3]

        The following flags can be used to define macros in your code.
        v=0|1 [0]         verbose                       VERBOSE 
        warn=0|1 [1]      warnings                      WARNINGS
        dbg=0|1 [0]                                     DSMDBG
        nruns= [1]                                      NRUNS

        D=[comma-sep list of host-only macros, without values]

        The options marked with * can be set as DSMBuilder.OPTION=VALUE in the SCons script
        The macros controlled by the other options are listed on the right
        The DSM Fortran API directory can be accessed via DSMBuilder.DSM_FORTRAN_API_DIR

        The environment variables listed in "$DSM_FORTRAN_API_DIR/dsm_env.sh" must be set.

        """         

    nth=getOpt('nth','#threads','1')
    loop_order=getOpt('order','loop order','1')
    nruns=getOpt('nruns','number of runs','1')

    env.Append(KERNEL_OPTS=['-DNRUNS='+nruns])
    env.Append(KERNEL_OPTS=['-DNTH='+nth])
    env.Append(KERNEL_OPTS=['-DLOOP_ORDER='+loop_order])
    ref=getOpt('ref','Reference','1')
    refflag='-DREF'
    verbose=getOpt('v','Verbose','0')
    vflag='-DVERBOSE'
    if verbose=='0':
        vflag=''

    warnings=getOpt('warn','Warnings','1')
    wflag='-Wall'
    if warnings=='0':
        vflag=''

    dbg=getOpt('dbg','Debug','0')    
    dbgflag='-g'
    dbgmacro='-DDSMDBG=1'
    if dbg=='0':
        dbgmacro=''
        dbgflag=''

    useF=getOpt('F','Functional',0)
    useFC='0'
    if 'useF' in env:
        if env['useF']==1:
            useF='1'
    if 'useFC' in env:
        if env['useFC']==1:
            useFC='1'
    useDyn=getOpt('dyn','Dynamic Library',0)
    if 'useDyn' in env:
        if env['useDyn']==1:
            useDyn='1'
            
    mcModel=getOpt('mcm','GCC Code Model Flag','s')
            
    optim=getOpt('O','Optimisation','3')
    optflag='-O'+optim
        
    defs=getOpt('D','Defines',None)
    defflags=[]
    if defs!=None:
        deflist=defs.split(',')
        defflags=map (lambda s: '-D'+s, deflist)   

    
    DEVFLAGS= env['KERNEL_OPTS']+ ['-DNRUNS='+nruns,'-DREF='+ref,vflag]+defflags
    if subprocess.check_output(['uname',''], shell=True).strip()  == "Darwin":
        DEVFLAGS+=['-DOSX']    
    env.Append(CXXFLAGS = ['-std=c++11',wflag,dbgflag,dbgmacro,optflag]+DEVFLAGS) 
    env.Append(CFLAGS = [wflag,dbgflag,optflag]+DEVFLAGS)     
    env['MACROS'] = DEVFLAGS
    #env.Append(CXXFLAGS = ['-mcmodel=large']

    env.Help(help)
    env.Append(CPPPATH=[DSM_FORTRAN_API_DIR+'/src'])
    env.Append(LIBS=['argo','argobackend-mpi'])

    if useF=='1':
        if 'FORTRAN_COMPILER' in os.environ:
            env['FORTRAN']=os.environ['FORTRAN_COMPILER']
            env['F95']=os.environ['FORTRAN_COMPILER']
        if 'FC' in os.environ:
            env['FORTRAN']=os.environ['FC']
            env['F95']=os.environ['FC']
        if ('GFORTRAN' in os.environ and env['FORTRAN'] == os.environ['GFORTRAN']) :
            env['FORTRANFLAGS']=env['CFLAGS']
            if OSX==1:
                env['LINKFLAGS']=['-Wl,-stack_size,0x40000000'] # Give OS X 1G stack
            env.Append(FORTRANFLAGS=['-Wno-aliasing','-Wno-unused','-Wno-unused-dummy-argument','-cpp','-m64','-ffree-form','-ffree-line-length-0','-fconvert=big-endian'])
            #env.Append(FORTRANFLAGS=['-mcmodel=large'])
#env['F95FLAGS']=['-Wno-aliasing','-Wno-unused','-Wno-unused-dummy-argument','-cpp','-m64','-mcmodel=medium','-ffree-form','-ffree-line-length-0','-fconvert=big-endian']
            env['F95FLAGS']=['-Wno-aliasing','-Wno-unused','-Wno-unused-dummy-argument','-cpp','-m64','-ffree-form','-ffree-line-length-0','-fconvert=big-endian']
            env.Append(F95FLAGS=env['CFLAGS'])
        else :
            env['CFLAGS'].pop(0)
            env['CFLAGS'].pop(0)
            env['CFLAGS'].pop(0)
            env['FORTRANFLAGS']=env['CFLAGS']
            if ('PGFORTRAN' in os.environ and env['FORTRAN'] == os.environ['PGFORTRAN']) : 
                env.Append(FORTRANFLAGS=['-cpp','-m64','-fast','-Mfree','-Mipa=fast'])
                env['F95FLAGS']=env['FORTRANFLAGS']
            else:
                print('Unknown compiler, no options specified.')
        if buildLib=='1':                
            if useDyn=='1':
                flib = env.SharedLibrary('dsmc', [DSM_FORTRAN_API_DIR+'/src/dsm-C.cc'])
            else:    
                flib = env.Library('dsmc', [DSM_FORTRAN_API_DIR+'/src/dsm-C.cc'])
            env.Depends('dsmAPI.o','dsmTypes.o')
            fflibt = env.Object('dsmTypes.o',DSM_FORTRAN_API_DIR+'/src/dsmTypes.f95')
            ffliba = env.Object('dsmAPI.o',DSM_FORTRAN_API_DIR+'/src/dsmAPI.f95')
            # env.Append(FORTRANMODDIR=[DSM_FORTRAN_API_DIR+'/lib/'])
            fflib = env.Library('dsmf', [fflibt,ffliba])
        else:
            flib = None
            fflib = None
            fflibt = None
            ffliba = None
    if useFC=='1':            
        if useDyn=='1':
            flib = env.SharedLibrary('dsmc', [DSM_FORTRAN_API_DIR+'/src/dsm-C.cc'])
        else:    
            flib = env.Library('dsmc', [DSM_FORTRAN_API_DIR+'/src/dsm-C.cc'])
        if 'MPICXX' in os.environ:
            env.Replace(CXX=os.environ['MPICXX'])
        else:
            env.Replace(CXX='mpicxx')
        if 'MPICC' in os.environ:
            env.Replace(CC=os.environ['MPICC'])
        else:
            env.Replace(CC='mpicc')
    env.Append(LIBPATH=['.',installPrefix+'/lib/'])
    if buildLib=='0':
        env.Append(FORTRANMODDIR=[installPrefix+'/lib/'])
    # env.Append(FORTRANPATH=['.',DSM_FORTRAN_API_DIR+'/src/'])
    if useF=='1':
        env.Append(LIBS=['stdc++'])
        if buildLib == '1':
            env.Install(installPrefix+'/lib', flib)
            env.Alias('install',installPrefix+'/lib', flib)
            env.Install(installPrefix+'/lib', fflib)
            env.Alias('installf',installPrefix+'/lib', fflib)
            env.Install(installPrefix+'/lib',['dsmapi.mod','dsmtypes.mod'])
        env.Prepend(LIBS=['dsmc'])
        env.Prepend(LIBS=['dsmf'])        

    elif useFC=='1':            
        env.Install(installPrefix+'/lib', flib)
        env.Alias('installf',installPrefix+'/lib', flib)
        env.Prepend(LIBS=['argobackend-mpi'])
        env.Prepend(LIBS=['argo'])
        env.Prepend(LIBS=['dsmc'])            

    else:
        env.Prepend(LIBS=['argobackend-mpi'])
        env.Prepend(LIBS=['argo'])

    return env

def build(appname,sources):
    global opts, DSM_FORTRAN_API_DIR
    env = initDSM()
    env.Program(appname,sources)

def buildF(env,appname,sources):
    initDSM(env)
    global opts, DSM_FORTRAN_API_DIR
    if useF=='1':
        env.Program(appname,sources)
    else:
        env.Program(appname,sources)

def buildLibs():
    envF=Environment(useF=1,buildLib=1)
    envF=initDSM(envF)