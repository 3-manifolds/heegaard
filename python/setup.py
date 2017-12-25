long_description =  """\
A Python package for using Berge's Heegaard to test if words can be
geometrically realized on the boundary of a handlebody.  
"""

import re, sys, subprocess, os, shutil, glob, sysconfig
from setuptools import setup, Command
from setuptools.command.build_py import build_py

def check_call(args):
    try:
        subprocess.check_call(args)
    except subprocess.CalledProcessError:
        executable = args[0]
        command = [a for a in args if not a.startswith('-')][-1]
        raise RuntimeError(command + ' failed for ' + executable)


class Clean(Command):
    """
    Removes the usual build/dist/egg-info directories.
    """
    user_options = []
    def initialize_options(self):
        pass 
    def finalize_options(self):
        pass
    def run(self):
        for dir in ['build', 'dist'] + glob.glob('*.egg-info'):
            if os.path.exists(dir):
                shutil.rmtree(dir)

class BuildPy(build_py):
    """
    Rebuilds the sqlite database files if needed.
    """
    def initialize_options(self):
        build_py.initialize_options(self)
        os.chdir('../src')        
        check_call(['make', 'heegaard_is_realizable'])
        os.chdir('../python')

class PipInstall(Command):
    user_options = [('no-build-wheel', 'n', 'assume wheel has already been built')]
    def initialize_options(self):
        self.no_build_wheel = False
    def finalize_options(self):
        pass
    def run(self): 
        python = sys.executable
        check_call([python, 'setup.py', 'build'])
        if not self.no_build_wheel:
            check_call([python, 'setup.py', 'bdist_wheel', '--universal'])
        egginfo = 'heegaard.egg-info'
        if os.path.exists(egginfo):
            shutil.rmtree(egginfo)
        wheels = glob.glob('dist' + os.sep + '*.whl')
        new_wheel = max(wheels, key=os.path.getmtime)            
        check_call([python, '-m', 'pip', 'install', '--upgrade',
                    '--upgrade-strategy', 'only-if-needed',
                    new_wheel])

class Test(Command):
    user_options = []
    def initialize_options(self):
        pass 
    def finalize_options(self):
        pass
    def run(self):
        build_lib_dir = os.path.join('build', 'lib')
        sys.path.insert(0, build_lib_dir)
        from heegaard.test import run_tests
        sys.exit(run_tests())

    
# Get version number from module
version = re.search("__version__ = '(.*)'",
                    open('src/__init__.py').read()).group(1)

setup(
    name = 'heegaard',
    version = version,
    description = 'checks realizability as a Heegaard diagram',
    long_description = long_description,
    url = 'https://bitbucket.org/t3m/heegaard',
    author = 'John Berge, Marc Culler, and Nathan M. Dunfield',
    author_email = 'snappy-help@computop.org',
    license='GPLv2+',
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
        'Operating System :: OS Independent',
        'Programming Language :: C',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Mathematics',
        ],

    packages = ['heegaard', 'heegaard/bin'],
    package_dir = {'heegaard':'src',
                   'heegaard/bin':'src/bin'},
    package_data = {'heegaard/bin': ['heegaard_is_realizable']},
    ext_modules = [],
    zip_safe = False,
    cmdclass = {'build_py': BuildPy,
                'clean': Clean,
                'pip_install':PipInstall,
                'test':Test
    },
)

