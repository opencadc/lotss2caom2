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


@patch('lotss2caom2.lotss_execute.LOTSSHierarchyStrategyContext._retrieve_provenance_metadata')
@patch('lotss2caom2.preview_augmentation.http_get')
@patch('lotss2caom2.lotss_execute.http_get')
@patch('lotss2caom2.clients.ASTRONClientCollection')
def test_preview_augmentation(clients_mock, http_get_mock, preview_get_mock, provenance_mock, test_config, tmp_path, test_name):
    import logging
    test_config.change_working_directory(tmp_path)
    clients_mock.py_vo_tap_client.search.side_effect = helpers._search_mosaic_id_mock

    def _endpoint_mock(url):
        result = type('response', (), {})()
        result.close = lambda: None
        with open(f'{test_name}/obs.xml') as f:
            result.text = f.read()
        return result

    clients_mock.https_session.get.side_effect = _endpoint_mock

    def _http_get_tar_mock(url, fqn, ignore_timeout):
        import logging
        logging.error(fqn)
        if url.endswith('preview=true'):
            shutil.copy(f'{test_name}/preview.jpg', fqn)
        else:
            shutil.copy(f'{test_name}/fits_headers.tar', fqn)

    http_get_mock.side_effect = _http_get_tar_mock
    preview_get_mock.side_effect = _http_get_tar_mock

    observation = read_obs_from_file(f'{test_name}/{basename(test_name)}_dr2.expected.xml')
    artifact_keys = get_all_artifact_keys(observation)
    if observation.observation_id == 'P124+62_dr2':
        assert len(artifact_keys) == 250, f'pre-condition artifact count {len(artifact_keys)}'
        provenance_mock.side_effect = helpers._get_db_query_mock_P164
    else:
        if observation.observation_id == 'P005+21_dr2':
            assert len(artifact_keys) == 493, f'pre-condition artifact count {len(artifact_keys)}'
        else:
            assert len(artifact_keys) == 7, f'pre-condition artifact count {len(artifact_keys)}'
        provenance_mock.return_value = []
    expander = lotss_execute.LOTSSHierarchyStrategyContext(clients_mock, test_config)
    expander.expand(test_name)
    test_config.working_directory = test_name
    for hierarchy in expander.hierarchies.values():
        hierarchy._working_directory = tmp_path
        hierarchy._preview_uri = f'https://vo.astron.nl/getproduct/LoTSS-DR2/{basename(test_name)}?preview=true'
        kwargs = {
            'hierarchy': hierarchy,
            'config': test_config,
        }
        observation = preview_augmentation.visit(observation, **kwargs)

    artifact_keys = get_all_artifact_keys(observation)
    if observation.observation_id == 'P124+62_dr2':
        assert len(artifact_keys) == 252, f'wrong number of artifacts {len(artifact_keys)}'
    else:
        if observation.observation_id == 'P005+21_dr2':
            assert len(artifact_keys) == 495, f'pre-condition artifact count {len(artifact_keys)}'
        else:
            assert len(artifact_keys) == 9, f'wrong number of artifacts {len(artifact_keys)}'
