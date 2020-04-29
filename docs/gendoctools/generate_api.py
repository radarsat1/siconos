"""Main driver to generate rst api/files for sphinx
in order to generate Siconos documentation.

Siconos is a program dedicated to modeling, simulation and control
 of non smooth dynamical systems.

 Copyright 2020 INRIA.

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
"""

from pathlib import Path
from gendoctools import cpp2rst
from gendoctools import python2rst


siconos_components = {'externals': 0, 'numerics': 1, 'kernel': 2, 'control': 3,
                      'mechanics': 4, 'mechanisms': 5, 'io': 6}


def build_rst_api(sphinx_directory):
    """Write api files for python and C++

    Parse existing rst files (one for each class,
    + those for functions)
    and collect them into cpp_api.rst and python_api.rst
    in sphinx/reference directory.

    This function is supposed to be called by a target
    generated with cmake (make rst_api).
    """

    # Place where files will be written
    outputdir = Path(sphinx_directory, 'reference')

    # Top file comments
    remark = '\n*If a file or a class you know does not appear in this page, '
    remark += 'it means it has not been (properly) documented or is not '
    remark += 'available in the high level API. Please contact us '
    remark += 'if you think it is an error.*\n\n'
    header = '.. contents::\n    :local:\n\n'
    introduction = 'Below you will find links to documentation'
    introduction += 'for all classes and files in Siconos, '
    introduction += 'sorted  by component.\n\n'
    common_header = remark + introduction + header
    # Create main rst file for C/C++ API
    cpp2rst.build_cpp_api_main(outputdir, common_header, siconos_components)
    # Create main rst file for Python API
    python2rst.build_python_api_main(outputdir, common_header,
                                     siconos_components)


def find_doxygen_diagrams(doxygen_path, output_directory):
    """Python utility to create sphinx rst file from png generated
    by doxygen (class diagrams ...)

    - scan doxygen ouput (html) path
    - create an image entry for each class_inh*.png file in a rst file

    Parameters
    ----------
    doxygen_path : Path()
       directory (full path) which contains png files generated by doxygen
    output_path : string
       directory (fullpath) of the rst output

    This function is supposed to be called  by a target created with
    cmake (make doxypng2sphinx)
    """

    # Scan doxygen output path and create a list with
    # files matching requirements
    class_diagram_match = 'inherit_graph*.png'
    doxygen_path = Path(doxygen_path).resolve()
    ref = '.. _api_class_diagrams:\n\n'
    header = 'C++ Class diagrams'
    header += '\n' + len(header) * '=' + '\n\n'
    header = ref + header
    files = [f for f in doxygen_path.glob(class_diagram_match)]
    realfiles = [Path('../doxygen', f.name) for f in files]
    realfiles = [Path('/', f) for f in realfiles]
    outputfile = Path(output_directory, 'class_diagrams.rst')
    with open(outputfile, 'w') as file:
        file.writelines(header)
        #params = [':height: 190 px', ':class: gallery']
        params = [':class: gallery']
        img_prefix = '.. image:: '

        for f in realfiles:
            line = img_prefix + f.as_posix()
            for p in params:
                line += '\n    ' + p
            line += '\n\n'
            file.writelines(line)
