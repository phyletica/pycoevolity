from setuptools import setup

from sumcoevolity import __version__

setup(
        name = "sumcoevolity",
        version = __version__,
        description = "Package for summarizing output of ecoevolity",
        author = "Jamie Oaks",
        author_email = "joaks1@gmail.com",
        license = "GPL",
        packages = ["sumcoevolity"],
        include_package_data = True,
        zip_safe = False,
        test_suite = "sumcoevolity.test.get_unittest_suite",
        # test_suite = "sumcoevolity.test",
        # test_loader = "unittest:TestLoader",
        install_requires = [
            # 'matplotlib'
        ],
        entry_points = {
            'console_scripts': [
                'loci2nex = sumcoevolity.cli.loci2nex:main',
            ],
        },
)
