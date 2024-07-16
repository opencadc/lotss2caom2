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

import logging
import traceback
from collections import deque

from mock import ANY, call, patch

from caom2pipe.manage_composable import read_obs_from_file, TaskType
from lotss2caom2 import composable, lotss_execute


@patch('lotss2caom2.clients.ASTRONClientCollection')
@patch('caom2pipe.execute_composable.OrganizeWithContext.do_one')
def test_run(run_mock, clients_mock, change_test_dir, test_config, tmp_path):
    run_mock.return_value = (0, None)
    test_f_id = 'test_file_id'
    test_f_name = f'{test_f_id}.fits'
    test_config.change_working_directory(tmp_path.as_posix())
    test_config.data_sources = ['pyvo_data_source:test']
    test_config.proxy_file_name = 'test_proxy.fqn'
    test_config.task_types = [TaskType.INGEST]
    test_config.write_to_file(test_config)

    with open(test_config.proxy_fqn, 'w') as f:
        f.write('test content')
    with open(test_config.work_fqn, 'w') as f:
        f.write(test_f_name)

    try:
        # execution
        test_result = composable._run()
    except Exception as e:
        logging.error(e)
        logging.error(traceback.format_exc())
        assert False, e

    assert test_result == 0, 'wrong return value'
    assert run_mock.called, 'should have been called'
    args, _ = run_mock.call_args
    assert args[0] == test_f_name, test_f_name


@patch('lotss2caom2.lotss_execute.ASTRONPyVODataSource')
@patch('lotss2caom2.clients.ASTRONClientCollection')
@patch('caom2pipe.execute_composable.OrganizeWithContext.do_one')
def test_run_remote_do_one_mock(run_mock, clients_mock, data_source_mock, change_test_dir, test_config, tmp_path):
    get_work_list = deque()
    test_obs_id = 'P000+23'
    get_work_list.append(test_obs_id)
    data_source_mock.return_value.get_work.return_value = get_work_list
    run_mock.return_value = (0, None)
    test_config.change_working_directory(tmp_path.as_posix())
    test_config.proxy_file_name = 'test_proxy.fqn'
    test_config.data_sources = ['pyvo_data_source:test']
    test_config.task_types = [TaskType.INGEST]
    test_config.write_to_file(test_config)

    with open(test_config.proxy_fqn, 'w') as f:
        f.write('test content')

    try:
        # execution
        test_result = composable._run_remote()
    except Exception as e:
        assert False, e

    assert test_result == 0, 'wrong return value'
    assert run_mock.called, 'should have been called'
    args, _ = run_mock.call_args
    assert args[0] == test_obs_id, test_obs_id


@patch('lotss2caom2.lotss_execute.LOTSSHierarchyStrategyContext.unset')
@patch('lotss2caom2.lotss_execute.LOTSSPreview.visit')
@patch('lotss2caom2.lotss_execute.LoTSSFits2caom2Visitor.visit')
@patch('lotss2caom2.lotss_execute.LOTSSHierarchyStrategyContext._expand')
@patch('lotss2caom2.lotss_execute.ASTRONClientCollection')
def test_run_client_mock(
    clients_mock,
    context_mock,
    visit_mock,
    preview_mock,
    unset_mock,
    test_data_dir,
    change_test_dir,
    test_config,
    tmp_path,
):
    # tests the circumstances where the entry driving the work is not the same as the observation ID for the tests
    test_obs_id = 'P000+23_dr2'
    test_uri = f'{test_config.collection}/{test_config.scheme}/{test_obs_id}/mosaic.fits'
    test_config.change_working_directory(tmp_path.as_posix())
    test_config.proxy_file_name = 'test_proxy.fqn'
    test_config.data_sources = ['pyvo_data_source:test']
    test_config.task_types = [TaskType.INGEST]
    test_config.logging_level = 'INFO'
    test_config.write_to_file(test_config)

    context_mock.return_value = {
        test_uri : lotss_execute.LOTSSHierarchyStrategy(test_uri, test_obs_id.replace("_dr2", "")),
    }
    clients_mock.return_value.metadata_client.read.return_value = None

    test_observation = read_obs_from_file(
        f'{test_data_dir}/{test_obs_id.replace("_dr2", "")}/{test_obs_id}.expected.xml'
    )
    visit_mock.return_value = test_observation
    preview_mock.return_value = test_observation

    with open(test_config.proxy_fqn, 'w') as f:
        f.write('test content')
    with open(test_config.work_fqn, 'w') as f:
        f.write(test_obs_id.replace("_dr2", ""))

    # execution
    test_result = composable._run()

    # post conditions
    assert test_result == 0, 'wrong return value'
    clients_mock.return_value.metadata_client.read.assert_has_calls(
        [call(test_config.collection, test_obs_id)],
    ), 'metadata read'
    clients_mock.return_value.metadata_client.create.assert_has_calls([call(ANY)]), 'metadata read'
    clients_mock.return_value.py_vo_tap_client.search.assert_has_calls([]), 'pyvo search'
    clients_mock.return_value.data_client.put.assert_has_calls([]), 'data put'
    clients_mock.return_value.data_client.get_head.assert_has_calls([]), 'data get_head'
