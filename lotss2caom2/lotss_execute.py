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

import logging
import tarfile
import traceback

from bs4 import BeautifulSoup
from datetime import datetime, timezone
from os import listdir, unlink
from os.path import basename, dirname, exists
from urllib.parse import urlparse

from caom2 import ProductType
from caom2pipe.astro_composable import get_vo_table_session, make_headers_from_file
from caom2pipe.data_source_composable import TodoFileDataSource
from caom2pipe.execute_composable import OrganizeWithContext
from caom2pipe.manage_composable import Config, http_get, Observable2, query_endpoint_session
from caom2pipe.run_composable import set_logging, TodoRunner
from caom2pipe.strategy_composable import HierarchyStrategy, HierarchyStrategyContext
from lotss2caom2.clients import ASTRONClientCollection
from lotss2caom2.data_source import ASTRONPyVODataSource
from lotss2caom2 import fits2caom2_augmentation
from lotss2caom2 import preview_augmentation


META_VISITORS = [fits2caom2_augmentation, preview_augmentation]
DATA_VISITORS = []


class LOTSSHierarchyStrategy(HierarchyStrategy):
    """
    The unit of work for a StorageName is the Observation ID. The file names are all found in the
    MetadataReader specialization based on that Observation ID.

    The destination URIs are set in the MetadataReader, and the file_uri is set to make preview generation work.

    Naming rules:
    - support mixed-case file name storage, and mixed-case obs id values
    - support uncompressed files in storage

    observationID example: P124+62
    """

    LOTSS_NAME_PATTERN = '*'

    def __init__(self, entry, mosaic_id):
        self._mosaic_id = mosaic_id
        self._uri_prefix = f'{HierarchyStrategy.scheme}:{HierarchyStrategy.collection}/{self._mosaic_id}'
        super().__init__(entry=entry, source_names=[entry], metadata=None, file_info=None)
        self._mosaic_metadata = None
        self._preview_uri = None

    def get_artifact_product_type(self, value):
        result = ProductType.SCIENCE
        if 'rms' in value:
            result = ProductType.NOISE
        elif 'pybdsmmask' in value or 'resid' in value:
            result = ProductType.AUXILIARY
        elif 'weights' in value:
            result = ProductType.WEIGHT
        elif '.jpg' in value:
            result = ProductType.PREVIEW
        return result

    @property
    def file_uri(self):
        return f'{self._uri_prefix}/{self._file_name}'

    @property
    def mosaic_id(self):
        return self._mosaic_id

    @property
    def mosaic_metadata(self):
        return self._mosaic_metadata

    def set_destination_uris(self):
        for entry in self._source_names:
            self._destination_uris.append(f'{self._uri_prefix}/{basename(entry)}')

    def set_product_id(self):
        if 'low-' in self._file_id:
            self._product_id = f'{self._mosaic_id}_mosaic_low'
        else:
            self._product_id = f'{self._mosaic_id}_mosaic'

    def set_obs_id(self):
        self._obs_id = f'{self._mosaic_id}_dr2'


class LOTSSRawHierarchyStrategy(HierarchyStrategy):
    """
    The unit of work for a StorageName is the Observation ID. The file names are all found in the
    MetadataReader specialization based on that Observation ID.

    The destination URIs are set in the MetadataReader, and the file_uri is set to make preview generation work.

    Naming rules:
    - support mixed-case file name storage, and mixed-case obs id values
    - support uncompressed files in storage

    observationID example: P124+62
    """

    LOTSS_NAME_PATTERN = '*'

    def __init__(self, mosaic_id, product_id, file_name):
        super().__init__(entry=file_name, source_names=[file_name], metadata=None, file_info=None)
        self._obs_id = f'{mosaic_id}_dr2'
        self._product_id = product_id
        self._logger.debug(self)

    def get_artifact_product_type(self, value):
        return ProductType.SCIENCE

    def set_obs_id(self, **kwargs):
        pass

    def set_product_id(self, **kwargs):
        pass


class LOTSSHierarchyStrategyContext(HierarchyStrategyContext):
    """This class takes one execution unit, and does the work, usually external to CADC, to make the HierarchyStratgy
    instances for it."""

    def __init__(self, clients, config):
        super().__init__(clients, config)
        self._service = clients.py_vo_tap_client
        self._session = clients.https_session
        self._mosaic_id = None
        self._mosaic_uri = None
        self._mosaic_metadata = None
        # the page that lists all the files for a MOSAIC
        self._related_products_uri = None
        # the Preview URL
        self._preview_uri = None
        # the URL that links to the tar of the FITS headers from the related_products URI
        self._headers_uri = None
        # the progenitor links to the raw data
        self._provenance_uris = []
        self._raw_table = {}
        self._http_get_timeout = config.http_get_timeout

    @property
    def mosaic_uri(self):
        return self._mosaic_uri

    def _get_mosaic_tap_metadata(self, mosaic_id):
        """Get all the metadata from the TAP query for a mosaic id

        Field names at time of writing:
        ('accref', 'owner', 'embargo', 'mime', 'accsize', 'centeralpha', 'centerdelta', 'imagetitle', 'instid',
         'dateobs', 'naxes', 'pixelsize', 'pixelscale', 'refframe', 'wcs_equinox', 'wcs_projection', 'wcs_refpixel',
         'wcs_refvalues', 'wcs_cdmatrix', 'bandpassid', 'bandpassunit', 'bandpassrefval', 'bandpasshi',
         'bandpasslo', 'pixflags', 'coverage', 'mosaic_id', 'related_products', 'lofar_obsids', 'data_pid',
         'adler32')
        """
        self._logger.debug(f'Begin _get_mosaic_tap_metadata for {mosaic_id}')
        results = self._service.search(f"SELECT * FROM lotss_dr2.mosaics WHERE mosaic_id='{mosaic_id}'")
        if len(results) == 1:
            self._mosaic_id = mosaic_id
            self._mosaic_uri = results[0]['accref']
            self._mosaic_metadata = results[0]
            strategy = LOTSSHierarchyStrategy(f'{self._mosaic_uri}/mosaic-blanked.fits', self._mosaic_id)
            strategy.metadata = [results[0]]
            strategy._mosaic_metadata = self._mosaic_metadata
            self._hierarchies[strategy.file_uri] = strategy
            self._related_products_uri = results[0]['related_products']
            self._logger.debug(f'Found mosaic metadata for {mosaic_id}')
        else:
            self._logger.warning(f'Wrong number of results {len(results)} when querying for mosaic {mosaic_id}.')
        self._logger.debug(f'End _get_mosaic_tap_metadata')

    def _get_related_products_metadata(self):
        """Go to the 'related_products' page and get all the metadata from that page."""
        self._logger.debug(f'Begin _get_related_products_metadata from {self._related_products_uri}')
        if self._related_products_uri:
            vo_table, error_message = get_vo_table_session(self._related_products_uri, self._session)
            if vo_table:
                for row in vo_table.array:
                    access_url = row['access_url']
                    # logging.error(access_url)
                    if 'preview' in access_url:
                        self._preview_uri = access_url
                        self._logger.debug(f'Found preview URL {access_url}')
                    elif 'fits' in access_url:
                        if 'headers.tar' in access_url:
                            self._headers_uri = access_url
                            self._logger.debug(f'Found Headers URL {access_url}')
                        else:
                            strategy = LOTSSHierarchyStrategy(access_url, self._mosaic_id)
                            strategy._mosaic_metadata = self._mosaic_metadata
                            self._hierarchies[strategy.file_uri] = strategy
                            # storage_name._destination_uris.append(access_url)
                            self._logger.debug(f'Found FITS URL {access_url}')
                    elif 'ObservationId=' in access_url:
                        self._provenance_uris.append(access_url)
                        self._logger.debug(f'Found Provenance at {access_url}')
            else:
                self._logger.warning(f'Encountered {error_message} when querying {self._related_products_uri}')
        self._logger.debug(f'End _get_related_products_metadata')

    def _get_headers_metadata(self):
        """Retrieve the file that has all the header metadata in it, and get all the metadata from that file,
        for each of mosaic'd files."""
        self._logger.debug(f'Begin _get_headers_metadata for {self._headers_uri}')
        if self._headers_uri:
            local_fqn = f'/tmp/{basename(self._headers_uri)}'
            if exists(local_fqn):
                unlink(local_fqn)
            http_get(self._headers_uri, local_fqn, self._http_get_timeout)
            if exists(local_fqn):
                with tarfile.open(local_fqn) as f:
                    f.extractall('/tmp', filter='data')
                for entry in listdir('/tmp/fits_headers'):
                    found_uri = None
                    if 'mosaic.rms.fits' in entry:
                        found_uri = f'{dirname(self._headers_uri)}/mosaic-rms.fits'
                    elif entry.endswith('.0.hdr'):
                        found_uri = f'{dirname(self._headers_uri)}/{basename(entry).replace(".0.hdr", "")}'
                    if found_uri:
                        temp_strategy = LOTSSHierarchyStrategy(found_uri, self._mosaic_id)
                        temp_strategy.metadata = make_headers_from_file(f'/tmp/fits_headers/{entry}')
                        temp_strategy._mosaic_metadata = self._mosaic_metadata
                        self._hierarchies[temp_strategy.file_uri] = temp_strategy
                        self._logger.info(f'Headers retrieved for {found_uri}')
                    else:
                        self._logger.warning(f'Unexpected file header {entry}')
        self._logger.debug('End _get_headers_metadata')

    def _get_provenance_metadata(self):
        """Retrieve available metadata for raw inputs."""
        if self._provenance_uris:
            for provenance_uri in self._provenance_uris:
                self._logger.debug(f'Begin _get_provenance_metadata from {provenance_uri}')
                self._raw_table = self._retrieve_provenance_metadata()
                sas_id = provenance_uri.split('=')[-1]
                for index, row in enumerate(self._raw_table):
                    if len(row) > 0:
                        if row.file_name:
                            strategy = LOTSSRawHierarchyStrategy(self._mosaic_id, sas_id, row.file_name)
                        else:
                            strategy = LOTSSRawHierarchyStrategy(self._mosaic_id, sas_id, str(index))
                        strategy.metadata = row
                        self._hierarchies[strategy.file_uri] = strategy
            self._logger.debug('End _get_provenance_metadata')

    def _retrieve_provenance_metadata(self):
        self._logger.debug(f'Begin _retrieve_provenance_metadata')
        response = None
        result = None
        try:
            response = query_endpoint_session(self._provenance_uri, self._session)
            response.raise_for_status()
            if response is None:
                self._logger.warning(f'No response from {self._provenance_uri}.')
            else:
                table = self._parse_html_string_for_table(response.content, 'result_table_AveragingPipeline')
                if table:
                    source_data_product = table[-1].get('Source DataProduct')
                    if source_data_product:
                        self._logger.info(f'Search {source_data_product} for raw metadata.')
                        response = query_endpoint_session(source_data_product, self._session)
                        response.raise_for_status()
                        if response is None:
                            self._logger.warning(f'No response from {source_data_product}')
                        else:
                            result = self._parse_html_string_for_table(
                                response.content, 'result_table_CorrelatedDataProduct'
                            )
                else:
                    self._logger.warning(f'Cannot find result_table_averagingPipeline in {self._provenance_uri}')
        finally:
            if response is not None:
                response.close()
        self._logger.debug(f'End _retrieve_provenance_metadata')
        return result

    def _expand(self, entry):
        self._logger.debug(f'Begin expand for {entry}')
        mosaic_id = basename(urlparse(entry).path)
        self._get_mosaic_tap_metadata(mosaic_id)
        self._get_related_products_metadata()
        self._get_headers_metadata()
        for hierarchy in self._hierarchies.values():
            hierarchy._preview_uri = self._preview_uri
        self._get_provenance_metadata()
        self._logger.debug(f'End expand with {len(self._hierarchies)} hierarchies.')
        return self._hierarchies

    def _parse_html_string_for_table(self, html_string, table_id):
        result = []
        soup = BeautifulSoup(html_string, 'html.parser')
        table = soup.find(id=table_id)
        if table:
            headers = [header.text.strip() for header in table.find_all('th')]
            # return [{headers[ii]: cell.text for ii, cell in enumerate(row.find_all('td'))} for row in table.find_all('tr')]
            for row in table.find_all('tr'):
                temp = {}
                for ii, cell in enumerate(row.find_all('td')):
                    if cell.css.select('a[href]'):
                        temp[headers[ii]] = cell.css.select('a[href]')[0].get('href')
                    else:
                        temp[headers[ii]] = cell.text
                result.append(temp)
        return result

    def set(self, entry):
        # this is here so that the CaomExecutor call doesn't fall over
        self._logger.debug(f'set for {entry}')


class StrategyTodoRunner(TodoRunner):

    def __init__(self, config, organizer, data_sources, observable):
        super().__init__(
            config,
            organizer,
            builder=None,
            data_sources=data_sources,
            metadata_reader=None,
            observable=observable,
            reporter=observable.reporter,
        )
        for data_source in data_sources:
            data_source.reporter = observable.reporter

    def _process_entry(self, data_source, entry, current_count):
        self._logger.debug(f'Begin _process_entry for {entry}.')
        start_s = datetime.now(tz=timezone.utc).timestamp()
        result, result_message = self._organizer.do_one(entry)
        try:
            data_source.clean_up(entry, result, current_count)
        except Exception as e:
            result_message = f'Data source cleanup failed with {e} for {entry}'
            self._logger.info(result_message)
            self._logger.debug(traceback.format_exc())
            result = -1
        if result == 0:
            self._reporter.capture_success_2(entry, start_s)
        else:
            self._reporter.capture_failure_2(entry, result_message)
        self._logger.debug(f'End _process_entry with result {result}.')
        return result

    def run(self):
        self._logger.debug('Begin run.')
        result = 0
        for data_source in self._data_sources:
            self._build_todo_list(data_source)
            # have the choose call here, so that retries don't change the set of tasks to be executed
            # self._organizer.choose()
            result |= self._run_todo_list(data_source, current_count=0)
        self._logger.debug('End run.')
        return result


def execute():
    logging.debug('Begin execute')
    # TODO - should I move the "get_executors" call into the Config constructor? they have to happen right after
    # each other anyway
    config = Config()
    config.get_executors()
    HierarchyStrategy.collection = config.collection
    HierarchyStrategy.preview_scheme = config.preview_scheme
    HierarchyStrategy.scheme = config.scheme
    HierarchyStrategy.data_source_extension = config.data_source_extensions
    set_logging(config)
    observable = Observable2(config)
    clients = ASTRONClientCollection(config)
    strategy_context = LOTSSHierarchyStrategyContext(clients, config)
    organizer = OrganizeWithContext(config, strategy_context, clients, observable)
    organizer.choose(META_VISITORS, DATA_VISITORS)
    data_source = TodoFileDataSource(config)
    # TODO - this would need to be consistent between data_sources
    runner = StrategyTodoRunner(
        config=config,
        organizer=organizer,
        data_sources=[data_source],
        observable=observable,
    )
    result = runner.run()
    result |= runner.run_retry()
    runner.report()
    logging.debug('End execute')
    return result


def remote_execute():
    logging.debug('Begin remote_execution')
    config = Config()
    config.get_executors()
    HierarchyStrategy.collection = config.collection
    HierarchyStrategy.preview_scheme = config.preview_scheme
    HierarchyStrategy.scheme = config.scheme
    HierarchyStrategy.data_source_extension = config.data_source_extensions
    set_logging(config)
    observable = Observable2(config)
    clients = ASTRONClientCollection(config)
    strategy_context = LOTSSHierarchyStrategyContext(clients, config)
    organizer = OrganizeWithContext(config, strategy_context, clients, observable)
    organizer.choose(META_VISITORS, DATA_VISITORS)
    data_source = ASTRONPyVODataSource(config, clients)
    runner = StrategyTodoRunner(
        config=config,
        organizer=organizer,
        data_sources=[data_source],
        observable=observable,
    )
    result = runner.run()
    result |= runner.run_retry()
    runner.report()
    logging.debug('End remote_execution')
    return result
