from setuptools import setup, find_packages

exec(open("pycoevolity/metadata.py").read())

setup(
        name = __project__,
        version = __version__,
        description = __description__,
        author = __author__,
        author_email = "joaks1@gmail.com",
        license = __license_short__,
        packages = find_packages(),
        include_package_data = True,
        zip_safe = False,
        test_suite = "pycoevolity.test.get_unittest_suite",
        # test_suite = "pycoevolity.test",
        # test_loader = "unittest:TestLoader",
        install_requires = [
            # 'matplotlib',
            'munkres >= 1.1.1; python_version > "3"',
            'munkres <= 1.0.12; python_version < "3"',
        ],
        entry_points = {
            'console_scripts': [
                'loci2nex = pycoevolity.cli.loci2alignment:main_nexus',
                'loci2phy = pycoevolity.cli.loci2alignment:main_phylip',
                'fastas2nex = pycoevolity.cli.fastas2alignment:main_nexus',
                'fastas2phy = pycoevolity.cli.fastas2alignment:main_phylip',
                'loci2dppmsbayes = pycoevolity.cli.loci2dppmsbayes:main',
                'fastas2fastas = pycoevolity.cli.fastas2fastas:main',
                'msb2nex = pycoevolity.cli.msbayes2alignments:main_nexus',
                'msb2phy = pycoevolity.cli.msbayes2alignments:main_phylip',
                'pyco-sumtimes = pycoevolity.cli.sumtimes:main',
                'pyco-sumevents = pycoevolity.cli.sumevents:main',
                'pyco-sumsizes = pycoevolity.cli.sumsizes:main',
                'pyco-sumchains = pycoevolity.cli.sumchains:main',
                'pyco-sumsims = pycoevolity.cli.sumsims:main',
                'pyco-dummy-data = pycoevolity.cli.dummy_data:main',
            ],
        },
)
