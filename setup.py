from setuptools import setup

from pycoevolity import __version__

setup(
        name = "pycoevolity",
        version = __version__,
        description = "Package for summarizing output of ecoevolity",
        author = "Jamie Oaks",
        author_email = "joaks1@gmail.com",
        license = "GPL",
        packages = ["pycoevolity"],
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
                'sumdivtimes = pycoevolity.cli.sumdivtimes:main',
                'sumpopsizes = pycoevolity.cli.sumpopsizes:main',
                'sumconvergence = pycoevolity.cli.sumconvergence:main',
            ],
        },
)
