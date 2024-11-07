# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2024.                            (c) 2024.
#  Government of Canada                 Gouvernement du Canada
#  National Research Council            Conseil national de recherches
#  Ottawa, Canada, K1A 0R6              Ottawa, Canada, K1A 0R6
#  All rights reserved                  Tous droits réservés
#
#  NRC disclaims any warranties,        Le CNRC dénie toute garantie
#  expressed, implied, or               énoncée, implicite ou légale,
#  statutory, of any kind with          de quelque nature que ce
#  respect to the software,             soit, concernant le logiciel,
#  including without limitation         y compris sans restriction
#  any warranty of merchantability      toute garantie de valeur
#  or fitness for a particular          marchande ou de pertinence
#  purpose. NRC shall not be            pour un usage particulier.
#  liable in any event for any          Le CNRC ne pourra en aucun cas
#  damages, whether direct or           être tenu responsable de tout
#  indirect, special or general,        dommage, direct ou indirect,
#  consequential or incidental,         particulier ou général,
#  arising from the use of the          accessoire ou fortuit, résultant
#  software.  Neither the name          de l'utilisation du logiciel. Ni
#  of the National Research             le nom du Conseil National de
#  Council of Canada nor the            Recherches du Canada ni les noms
#  names of its contributors may        de ses  participants ne peuvent
#  be used to endorse or promote        être utilisés pour approuver ou
#  products derived from this           promouvoir les produits dérivés
#  software without specific prior      de ce logiciel sans autorisation
#  written permission.                  préalable et particulière
#                                       par écrit.
#
#  This file is part of the             Ce fichier fait partie du projet
#  OpenCADC project.                    OpenCADC.
#
#  OpenCADC is free software:           OpenCADC est un logiciel libre ;
#  you can redistribute it and/or       vous pouvez le redistribuer ou le
#  modify it under the terms of         modifier suivant les termes de
#  the GNU Affero General Public        la “GNU Affero General Public
#  License as published by the          License” telle que publiée
#  Free Software Foundation,            par la Free Software Foundation
#  either version 3 of the              : soit la version 3 de cette
#  License, or (at your option)         licence, soit (à votre gré)
#  any later version.                   toute version ultérieure.
#
#  OpenCADC is distributed in the       OpenCADC est distribué
#  hope that it will be useful,         dans l’espoir qu’il vous
#  but WITHOUT ANY WARRANTY;            sera utile, mais SANS AUCUNE
#  without even the implied             GARANTIE : sans même la garantie
#  warranty of MERCHANTABILITY          implicite de COMMERCIALISABILITÉ
#  or FITNESS FOR A PARTICULAR          ni d’ADÉQUATION À UN OBJECTIF
#  PURPOSE.  See the GNU Affero         PARTICULIER. Consultez la Licence
#  General Public License for           Générale Publique GNU Affero
#  more details.                        pour plus de détails.
#
#  You should have received             Vous devriez avoir reçu une
#  a copy of the GNU Affero             copie de la Licence Générale
#  General Public License along         Publique GNU Affero avec
#  with OpenCADC.  If not, see          OpenCADC ; si ce n’est
#  <http://www.gnu.org/licenses/>.      pas le cas, consultez :
#                                       <http://www.gnu.org/licenses/>.
#
#  Revision: 4
#
# ***********************************************************************
#

from caom2pipe.manage_composable import read_obs_from_file, write_obs_to_file
from lotss2caom2.lotss_execute import LOTSSHierarchyStrategy
from lotss2caom2.position_bounds_augmentation import visit


def test_visit(test_data_dir, test_config, tmp_path):
    test_config.change_working_directory(tmp_path)
    test_mosaic_id = 'P002+18'
    test_dir = f'{test_data_dir}/{test_mosaic_id}'
    observation = read_obs_from_file(f'{test_dir}/{test_mosaic_id}_dr2.expected.xml')
    test_config.working_directory = tmp_path
    test_fqn = f'{test_dir}/mosaic-blanked.fits'
    hierarchy = LOTSSHierarchyStrategy(test_fqn, test_mosaic_id)
    test_artifact = observation.planes[hierarchy.product_id].artifacts[hierarchy.file_uri]
    # for part in test_artifact.parts.values():
    #     for chunk in part.chunks:
    #         assert chunk.position.axis.bounds is None, f'{part.name} precondition'
    kwargs = {
        'hierarchy': hierarchy,
        'working_directory': f'{test_data_dir}/{test_mosaic_id}/{test_mosaic_id}',
        'log_file_directory': None,
    }
    observation = visit(observation, **kwargs)
    assert observation is not None, 'expect a return value'
    for part in test_artifact.parts.values():
        for chunk in part.chunks:
            assert chunk.position.axis.bounds is not None, f'{part.name} postcondition'

    write_obs_to_file(observation, './fpf_end.xml')

