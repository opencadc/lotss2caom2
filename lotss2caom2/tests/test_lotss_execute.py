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
#  $Revision: 4 $
#
# ***********************************************************************
#

import shutil

from caom2pipe.caom_composable import get_all_artifact_keys
from caom2pipe.execute_composable import MetaVisitExpander, OrganizeWithContext, OrganizeWithHierarchy
from caom2pipe.execute_composable import MetaVisitHierarchies
from caom2pipe.manage_composable import Config, Observable2, TaskType
from lotss2caom2.lotss_execute import DATA_VISITORS, execute, LOTSSHierarchyStrategyContext, META_VISITORS, remote_execute
from lotss2caom2.lotss_execute import DATA_VISITORS_FASTER, META_VISITORS_FASTER
from mock import ANY, call, patch
from unittest import skip

import helpers


@patch('lotss2caom2.lotss_execute.LOTSSHierarchyStrategyContext._retrieve_provenance_metadata')
@patch('lotss2caom2.lotss_execute.http_get')
@patch('lotss2caom2.clients.ASTRONClientCollection')
def test_strategy(clients_mock, http_get_mock, session_mock, test_config, test_data_dir, tmp_path):
    # url2 - P000+23 = https://lta.lofar.eu/Lofar?project=ALL&product=all_observation_pipeline&mode=query_result_page_user&ObservationId=689778
    # url3 - P000+23 = https://lta.lofar.eu/Lofar?project=ALL&mode=query_result_page&product=CorrelatedDataProduct&pipeline_object_id=7CB88D2B6B6F3778E053164A17ACB062

    test_config.change_working_directory(tmp_path)
    clients_mock.py_vo_tap_client.search.side_effect = helpers._search_mosaic_id_mock

    def _endpoint_mock(url):
        assert (
            url == 'https://vo.astron.nl/lotss_dr2/q/dlmosaic/dlmeta?ID=ivo%3A//astron.nl/%7E%3FLoTSS-DR2/P000%2B23'
        ), f'wrong url {url}'
        result = type('response', (), {})()
        result.close = lambda: None
        with open(f'{test_data_dir}/P000+23/obs.xml') as f:
            result.text = f.read()
        return result

    clients_mock.https_session.get.side_effect = _endpoint_mock

    def _http_get_mock(url, fqn, ignore_timeout):
        assert url == 'https://lofar-webdav.grid.surfsara.nl:2881/P000+23/fits_headers.tar', f'wrong url {url}'
        assert fqn == '/tmp/fits_headers.tar', f'wrong url {fqn}'
        shutil.copy(f'{test_data_dir}/P000+23/fits_headers.tar', '/tmp')

    http_get_mock.side_effect = _http_get_mock

    session_mock.side_effect = helpers._get_db_query_mock_P164

    test_subject = LOTSSHierarchyStrategyContext(clients_mock, test_config)
    assert test_subject is not None, 'ctor'

    test_mosaic_id = 'P000+23'
    test_subject.expand(test_mosaic_id)

    assert len(test_subject.hierarchies) == 250, f'wrong header length {len(test_subject.hierarchies)}'
    assert 'astron:LOTSS/P000+23/mosaic.fits' not in test_subject.hierarchies.keys(), 'the renaming bug is back'

    test_uri_prefix = f'{test_config.scheme}:{test_config.collection}/{test_mosaic_id}/'
    for check in [test_subject.hierarchies.keys()]:
        assert f'{test_uri_prefix}low-mosaic-blanked.fits' in check, f'{check.__class__.__name__} low-mosaic-blanked'
        assert f'{test_uri_prefix}low-mosaic-weights.fits' in check, f'{check.__class__.__name__} low-mosaic-weights'
        assert f'{test_uri_prefix}mosaic-blanked.fits' in check, f'{check.__class__.__name__} mosaic-blanked'
        assert f'{test_uri_prefix}mosaic-weights.fits' in check, f'{check.__class__.__name__} mosaic-weights'
        assert f'{test_uri_prefix}mosaic.pybdsmmask.fits' in check, f'{check.__class__.__name__} mosaic.pybdsmmask'
        assert f'{test_uri_prefix}mosaic.resid.fits' in check, f'{check.__class__.__name__} mosaic.resid'
        assert f'{test_uri_prefix}mosaic-rms.fits' in check, f'{check.__class__.__name__} mosaic.rms'
        assert f'{test_config.scheme}:{test_config.collection}/L664568_SB484_uv.MS_614524cc.tar' in check, f'{check.__class__.__name__} MS'


@patch('lotss2caom2.lotss_execute.query_endpoint_session')
@patch('lotss2caom2.lotss_execute.http_get')
@patch('lotss2caom2.clients.ASTRONClientCollection')
def test_strategy_https(clients_mock, http_get_mock, session_mock, test_config, test_data_dir, tmp_path):
    test_config.change_working_directory(tmp_path)
    clients_mock.py_vo_tap_client.search.side_effect = helpers._search_mosaic_id_mock

    def _endpoint_mock(url):
        assert (
            url == 'https://vo.astron.nl/lotss_dr2/q/dlmosaic/dlmeta?ID=ivo%3A//astron.nl/%7E%3FLoTSS-DR2/P000%2B23'
        ), f'wrong url {url}'
        result = type('response', (), {})()
        result.close = lambda: None
        with open(f'{test_data_dir}/P000+23/obs.xml') as f:
            result.text = f.read()
        return result

    clients_mock.https_session.get.side_effect = _endpoint_mock

    def _http_get_mock(url, fqn, ignore_timeout):
        assert url == 'https://lofar-webdav.grid.surfsara.nl:2881/P000+23/fits_headers.tar', f'wrong url {url}'
        assert fqn == '/tmp/fits_headers.tar', f'wrong url {fqn}'
        shutil.copy(f'{test_data_dir}/P000+23/fits_headers.tar', '/tmp')

    http_get_mock.side_effect = _http_get_mock

    def _session_mock(url, ignore_session):
        # first if - P000+23 = https://lta.lofar.eu/Lofar?project=ALL&product=all_observation_pipeline&
        #   mode=query_result_page_user&ObservationId=689778
        # second if - P000+23 = https://lta.lofar.eu/Lofar?project=ALL&mode=query_result_page&
        #   product=CorrelatedDataProduct&pipeline_object_id=7CB88D2B6B6F3778E053164A17ACB062
        result = type('response', (), {})()
        result.close = lambda: None
        result.raise_for_status = lambda: None
        if url == 'https://lta.lofar.eu/Lofar?project=ALL&product=all_observation_pipeline&mode=query_result_page_user&ObservationId=689778':
            with open(f'{test_data_dir}/provenance/p000+23_averaging_pipeline.html') as f:
                result.content = f.read()
        elif url == 'https://lta.lofar.eu/Lofar?project=ALL&mode=query_result_page&product=CorrelatedDataProduct&pipeline_object_id=7CB88D2B6B6F3778E053164A17ACB062':
            with open(f'{test_data_dir}/provenance/p000+23_correlated_data_products.html') as f:
                result.content = f.read()
        else:
            assert False, url
        return result

    session_mock.side_effect = _session_mock

    test_subject = LOTSSHierarchyStrategyContext(clients_mock, test_config)
    assert test_subject is not None, 'ctor'

    test_mosaic_id = 'P000+23'
    test_subject.expand(test_mosaic_id)

    assert len(test_subject.hierarchies) == 250, f'wrong header length {len(test_subject.hierarchies)}'
    assert 'astron:LOTSS/P000+23/mosaic.fits' not in test_subject.hierarchies.keys(), 'the renaming bug is back'

    test_uri_prefix = f'{test_config.scheme}:{test_config.collection}/{test_mosaic_id}/'
    for check in [test_subject.hierarchies.keys()]:
        assert f'{test_uri_prefix}low-mosaic-blanked.fits' in check, f'{check.__class__.__name__} low-mosaic-blanked'
        assert f'{test_uri_prefix}low-mosaic-weights.fits' in check, f'{check.__class__.__name__} low-mosaic-weights'
        assert f'{test_uri_prefix}mosaic-blanked.fits' in check, f'{check.__class__.__name__} mosaic-blanked'
        assert f'{test_uri_prefix}mosaic-weights.fits' in check, f'{check.__class__.__name__} mosaic-weights'
        assert f'{test_uri_prefix}mosaic.pybdsmmask.fits' in check, f'{check.__class__.__name__} mosaic.pybdsmmask'
        assert f'{test_uri_prefix}mosaic.resid.fits' in check, f'{check.__class__.__name__} mosaic.resid'
        assert f'{test_uri_prefix}mosaic-rms.fits' in check, f'{check.__class__.__name__} mosaic.rms'
        assert f'{test_config.scheme}:{test_config.collection}/L689778_SB006_uv.MS_fb14a2ee.tar' in check, f'{check.__class__.__name__} MS'


@patch('lotss2caom2.fits2caom2_augmentation.visit')
@patch('lotss2caom2.preview_augmentation.visit')
@patch('lotss2caom2.lotss_execute.LOTSSHierarchyStrategyContext._retrieve_provenance_metadata')
@patch('lotss2caom2.lotss_execute.http_get')
@patch('lotss2caom2.clients.ASTRONClientCollection')
def test_organize_nominal(
    clients_mock,
    http_get_mock,
    session_mock,
    preview_mock,
    visit_mock,
    test_config,
    test_data_dir,
    tmp_path,
    change_test_dir,
):
    # one mosaic id, seven files
    test_config.change_working_directory(tmp_path)
    clients_mock.py_vo_tap_client.search.side_effect = helpers._search_mosaic_id_mock

    def _endpoint_mock(url):
        result = type('response', (), {})()
        result.close = lambda: None
        with open(f'{test_data_dir}/P000+23/obs.xml') as f:
            result.text = f.read()
        return result

    clients_mock.https_session.get.side_effect = _endpoint_mock

    def _http_get_mock(url, fqn, ignore_timeout):
        shutil.copy(f'{test_data_dir}/P000+23/fits_headers.tar', '/tmp')

    http_get_mock.side_effect = _http_get_mock

    session_mock.side_effect = helpers._get_db_query_mock
    clients_mock.metadata_client.read.side_effect = helpers._observation

    def _visit_mock(observation, **kwargs):
        import logging
        logging.error(f'yes, no??????? {observation}')
        return observation

    preview_mock.side_effect = _visit_mock
    visit_mock.side_effect = _visit_mock

    test_context = LOTSSHierarchyStrategyContext(clients_mock, test_config)

    test_config.task_types = [TaskType.INGEST]
    test_config.change_working_directory(tmp_path)
    test_observable = Observable2(test_config)
    test_subject = OrganizeWithContext(test_config, test_context, clients_mock, test_observable)
    test_subject.choose(META_VISITORS, DATA_VISITORS)
    assert len(test_subject._executors) == 1, len(test_subject._executors)
    assert isinstance(test_subject._executors[0], MetaVisitExpander), test_subject._executors
    for entry in [
        'file_name.fits',
        'cadc:COLLECTION/file_name.fits',
        'https://a/long/path/file_name.fits',
        'pawsey9989:possum/tiles',
        'observationID',
    ]:
        test_result_value, test_result_message = test_subject.do_one(entry)
        assert test_result_value == 0, f'expect success {entry} {test_result_message}'
        assert test_result_message is None, f'success means no message {entry} {test_result_message}'


@skip('not implemented')
def test_organize_failures():
    # set up each of the failure paths
    assert False


@skip('not properly mocked')
@patch('lotss2caom2.lotss_execute.ASTRONClientCollection')
def test_remote_execute_nominal(clients_mock, test_config, tmp_path, change_test_dir):
    test_config.change_working_directory(tmp_path)
    test_config.proxy_file_name = 'test_proxy.pem'
    test_config.data_sources = ['https://localhost']
    test_config.task_types = [TaskType.INGEST]
    test_config.logging_level = 'DEBUG'

    with open(test_config.proxy_fqn, 'w') as f:
        f.write('test content')

    Config.write_to_file(test_config)

    clients_mock.return_value.py_vo_tap_client.search.side_effect = [
        helpers._search_id_list_mock(None),
        helpers._search_mosaic_id_mock('P000'),
        helpers._search_mosaic_id_mock('P000'),
        helpers._search_mosaic_id_mock('P000'),
    ]
    clients_mock.return_value.https_session.get.side_effect = helpers._get_vo_mock

    test_result = remote_execute()
    assert test_result == 0, 'expect success'
    clients_mock.return_value.py_vo_tap_client.search.assert_has_calls(
        [
            call('SELECT mosaic_id from lotss_dr2.mosaics'),
            call("SELECT * FROM lotss_dr2.mosaics WHERE mosaic_id='P000+23'"),
            call("SELECT * FROM lotss_dr2.mosaics WHERE mosaic_id='P000+31'"),
            call("SELECT * FROM lotss_dr2.mosaics WHERE mosaic_id='P000+36'"),
        ],
    ), 'py_vo search'
    clients_mock.return_value.data_client.put.assert_has_calls(
        []
    ), f'data {clients_mock.return_value.data_client.put.calls}'
    clients_mock.return_value.metadata_client.create.assert_has_calls([]), 'metadata create'
    clients_mock.return_value.metadata_client.read.assert_has_calls([]), 'metadata read'
    assert False


@patch('lotss2caom2.preview_augmentation.http_get')
@patch('lotss2caom2.lotss_execute.query_endpoint_session')
@patch('lotss2caom2.lotss_execute.http_get')
@patch('lotss2caom2.lotss_execute.ASTRONClientCollection')
def test_execute_nominal(
    clients_mock,
    lotss_execute_http_get_mock,
    session_mock,
    preview_http_get_mock,
    test_config,
    tmp_path,
    change_test_dir,
):
    test_config.change_working_directory(tmp_path)
    test_config.proxy_file_name = 'test_proxy.pem'
    test_config.task_types = [TaskType.INGEST]

    with open(test_config.proxy_fqn, 'w') as f:
        f.write('test content')
    with open(test_config.work_fqn, 'w') as f:
        f.write('test_content')

    Config.write_to_file(test_config)

    clients_mock.return_value.py_vo_tap_client.search.side_effect = helpers._search_mosaic_id_mock
    clients_mock.return_value.metadata_client.read.side_effect = lambda collection, obs_id: None
    test_result = execute()
    assert test_result == 0, 'expect success'
    clients_mock.return_value.py_vo_tap_client.search.assert_called_with(
        "SELECT * FROM lotss_dr2.mosaics WHERE mosaic_id='test_content'",
    ), 'py_vo search'
    assert not clients_mock.return_value.data_client.read.called, 'data'
    clients_mock.return_value.metadata_client.read.assert_called_with('LOTSS', 'test_content_dr2'), 'metadata'
    clients_mock.return_value.metadata_client.create.assert_called_with(ANY), 'metadata create'


@patch('lotss2caom2.fits2caom2_augmentation_faster.visit')
@patch('lotss2caom2.preview_augmentation_faster.visit')
@patch('lotss2caom2.lotss_execute.LOTSSHierarchyStrategyContext._retrieve_provenance_metadata')
@patch('lotss2caom2.lotss_execute.http_get')
@patch('lotss2caom2.clients.ASTRONClientCollection')
def test_organize_faster_nominal(
    clients_mock,
    http_get_mock,
    session_mock,
    preview_mock,
    visit_mock,
    test_config,
    test_data_dir,
    tmp_path,
    change_test_dir,
):
    # one mosaic id, 493 files
    test_config.change_working_directory(tmp_path)
    clients_mock.py_vo_tap_client.search.side_effect = helpers._search_mosaic_id_mock

    def _endpoint_mock(url):
        result = type('response', (), {})()
        result.close = lambda: None
        with open(f'{test_data_dir}/P005+21/obs.xml') as f:
            result.text = f.read()
        return result

    clients_mock.https_session.get.side_effect = _endpoint_mock

    def _http_get_mock(url, fqn, ignore_timeout):
        shutil.copy(f'{test_data_dir}/P005+21/fits_headers.tar', '/tmp')

    http_get_mock.side_effect = _http_get_mock

    session_mock.side_effect = [helpers._get_db_query_mock_P005A(), helpers._get_db_query_mock_P005B()]
    clients_mock.metadata_client.read.side_effect = helpers._observation

    def _visit_mock(observation, **kwargs):
        import logging
        logging.error('visit mock!!!!!')
        return observation

    preview_mock.side_effect = _visit_mock
    visit_mock.side_effect = _visit_mock

    test_context = LOTSSHierarchyStrategyContext(clients_mock, test_config)

    test_config.task_types = [TaskType.INGEST]
    test_config.change_working_directory(tmp_path)
    test_observable = Observable2(test_config)
    test_subject = OrganizeWithHierarchy(test_config, test_context, clients_mock, test_observable)
    test_subject.choose(META_VISITORS_FASTER, DATA_VISITORS_FASTER)
    assert len(test_subject._executors) == 1, len(test_subject._executors)
    assert isinstance(test_subject._executors[0], MetaVisitHierarchies), test_subject._executors

    entry = 'P005+21'
    test_result_value, test_result_message = test_subject.do_one(entry)
    assert test_result_value == 0, f'expect success {entry} {test_result_message}'
    assert test_result_message is None, f'success means no message {entry} {test_result_message}'
    assert visit_mock.called, 'visit called'
    assert visit_mock.call_count == 1, 'visit call count'
    assert preview_mock.called, 'preview visit called'
    assert preview_mock.call_count == 1, 'preview visit call count'
