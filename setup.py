from setuptools import setup, find_packages

from pycoevolity import __version__

setup(
        name = "pycoevolity",
        version = __version__,
        description = "Package for summarizing output of ecoevolity",
        author = "Jamie Oaks",
        author_email = "joaks1@gmail.com",
        license = "GPL",
        packages = find_packages(),
        include_package_data = True,
        zip_safe = False,
        test_suite = "pycoevolity.test.get_unittest_suite",
        # test_suite = "pycoevolity.test",
        # test_loader = "unittest:TestLoader",
        install_requires = [
            # 'matplotlib'
        ],
        entry_points = {
            'console_scripts': [
                'loci2nex = pycoevolity.cli.loci2alignment:main_nexus',
                'loci2phy = pycoevolity.cli.loci2alignment:main_phylip',
                'fastas2nex = pycoevolity.cli.fastas2alignment:main_nexus',
                'fastas2phy = pycoevolity.cli.fastas2alignment:main_phylip',
                'loci2dppmsbayes = pycoevolity.cli.loci2dppmsbayes:main',
                'fastas2fastas = pycoevolity.cli.fastas2fastas:main',
                'pyco-sumtimes = pycoevolity.cli.sumtimes:main',
                'pyco-sumevents = pycoevolity.cli.sumevents:main',
                'pyco-sumsizes = pycoevolity.cli.sumsizes:main',
                'pyco-sumchains = pycoevolity.cli.sumchains:main',
                'pyco-sumsims = pycoevolity.cli.sumsims:main',
            ],
        },
)
