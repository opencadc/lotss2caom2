# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2023.                            (c) 2023.
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

import glob
import shutil
from os.path import basename
from caom2pipe.caom_composable import get_all_artifact_keys
from caom2pipe.manage_composable import read_obs_from_file
from lotss2caom2 import preview_augmentation
from lotss2caom2 import lotss_execute

from mock import Mock, patch
import helpers


def pytest_generate_tests(metafunc):
    obs_id_list = glob.glob(f'{metafunc.config.invocation_dir}/data/P*')
    metafunc.parametrize('test_name', obs_id_list)


@patch('lotss2caom2.preview_augmentation.http_get')
def test_preview_augmentation(preview_get_mock, test_config, tmp_path, test_name):
    test_config.change_working_directory(tmp_path)

    def _http_get_tar_mock(url, fqn, _):
        import logging
        logging.error(fqn)
        if url.endswith('preview=true'):
            shutil.copy(f'{test_name}/preview.jpg', fqn)
        else:
            shutil.copy(f'{test_name}/fits_headers.tar', fqn)

    preview_get_mock.side_effect = _http_get_tar_mock
    clients_mock = Mock()

    observation = read_obs_from_file(f'{test_name}/{basename(test_name)}_dr2.expected.xml')
    artifact_keys = get_all_artifact_keys(observation)
    assert len(artifact_keys) == 7, f'pre-condition artifact count {len(artifact_keys)}'
    hierarchies = {
        f'{test_config.scheme}:{test_config.collection}/{basename(test_name)}/mosaic.fits':
            lotss_execute.LOTSSHierarchyStrategy(entry='mosaic.fits', mosaic_id=f'{basename(test_name)}'),
        f'{test_config.scheme}:{test_config.collection}/{basename(test_name)}/low-mosaic-blanked.fits':
            lotss_execute.LOTSSHierarchyStrategy(entry='low-mosaic-blanked.fits', mosaic_id=f'{basename(test_name)}'),
    }
    expander = lotss_execute.LOTSSHierarchyStrategyContext(clients_mock, test_config)
    expander.expand(test_name)
    test_config.working_directory = test_name
    for hierarchy in hierarchies.values():
        hierarchy._working_directory = tmp_path
        hierarchy._preview_uri = f'https://vo.astron.nl/getproduct/LoTSS-DR2/{basename(test_name)}?preview=true'
        kwargs = {
            'hierarchy': hierarchy,
            'config': test_config,
        }
        observation = preview_augmentation.visit(observation, **kwargs)

    artifact_keys = get_all_artifact_keys(observation)
    assert len(artifact_keys) == 11, f'wrong number of artifacts {len(artifact_keys)}'
