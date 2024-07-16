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
#  $Revision: 4 $
#
# ***********************************************************************
#

from mock import patch

from lotss2caom2 import fits2caom2_augmentation
from caom2.diff import get_differences
from caom2pipe import manage_composable as mc
from lotss2caom2 import lotss_execute

import glob
import helpers
import os
import shutil

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_DATA_DIR = os.path.join(THIS_DIR, 'data')


def pytest_generate_tests(metafunc):
    obs_id_list = glob.glob(f'{TEST_DATA_DIR}/P*')
    # obs_id_list = glob.glob(f'{TEST_DATA_DIR}/P000+23')
    metafunc.parametrize('test_name', obs_id_list)


@patch('lotss2caom2.lotss_execute.LOTSSHierarchyStrategyContext._retrieve_provenance_metadata')
@patch('lotss2caom2.lotss_execute.query_endpoint_session')
@patch('lotss2caom2.lotss_execute.http_get')
@patch('lotss2caom2.clients.ASTRONClientCollection')
def test_main_app(
    clients_mock,
    http_get_mock,
    session_mock,
    retrieve_provenance_mock,
    test_name,
    test_data_dir,
    test_config,
    tmp_path,
):
    test_config.change_working_directory(tmp_path)
    clients_mock.py_vo_tap_client.search.side_effect = helpers._search_mosaic_id_mock

    def _endpoint_mock(url):
        result = type('response', (), {})()
        result.close = lambda: None
        with open(f'{test_name}/obs.xml') as f:
            result.text = f.read()
        return result

    clients_mock.https_session.get.side_effect = _endpoint_mock

    def _http_get_mock(url, fqn, ignore_timeout):
        assert fqn == '/tmp/fits_headers.tar', f'wrong url {fqn}'
        shutil.copy(f'{test_name}/fits_headers.tar', '/tmp')

    http_get_mock.side_effect = _http_get_mock

    def _session_mock(url, _):
        result = type('response', (), {})()
        result.close = lambda: None
        result.raise_for_status = lambda: None
        if url == 'https://lta.lofar.eu/Lofar?project=ALL&product=all_observation_pipeline&mode=query_result_page_user&ObservationId=689778':
            with open(f'{test_data_dir}/provenance/progenitor.html') as f:
                result.content = f.read()
        else:
            with open(f'{test_data_dir}/provenance/source_data_products.html') as f:
                result.content = f.read()
        return result
    session_mock.side_effect = _session_mock

    retrieve_provenance_mock.side_effect = helpers._get_db_query_mock

    expected_fqn = f'{test_name}/{os.path.basename(test_name)}_dr2.expected.xml'
    actual_fqn = expected_fqn.replace('expected', 'actual')
    if os.path.exists(actual_fqn):
        os.unlink(actual_fqn)

    observations = {}
    expander = lotss_execute.LOTSSHierarchyStrategyContext(clients_mock, test_config)
    expander.expand(test_name)
    for hierarchy in expander.hierarchies.values():
        kwargs = {
            'hierarchy': hierarchy,
            'config': test_config,
        }
        observation = observations.get(hierarchy.obs_id)
        visitor  = fits2caom2_augmentation.LoTSSFits2caom2Visitor()
        observation = visitor.visit(observation, **kwargs)
        observations[hierarchy.obs_id] = observation

    if len(observations) == 0:
        assert False, f'Did not create observation for {test_name}'
    else:
        for observation in observations.values():
            if os.path.exists(expected_fqn):
                expected = mc.read_obs_from_file(expected_fqn)
                compare_result = get_differences(expected, observation)
                if compare_result is not None:
                    mc.write_obs_to_file(observation, actual_fqn)
                    compare_text = '\n'.join([r for r in compare_result])
                    msg = f'Differences found in observation {expected.observation_id}\n' f'{compare_text}'
                    raise AssertionError(msg)
            else:
                mc.write_obs_to_file(observation, actual_fqn)
                assert False, f'nothing to compare to for {test_name}, missing {expected_fqn}'
    # assert False  # cause I want to see logging messages
