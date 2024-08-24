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

import logging
import os
import shutil

from caom2 import ReleaseType, ProductType
from caom2pipe.manage_composable import CadcException, get_artifact_metadata, http_get, PreviewMeta, PreviewVisitor


def search_for_file(strategy, file_name, working_directory):
    temp_fqn = os.path.join(working_directory, file_name)
    temp_obs_fqn = os.path.join(os.path.join(working_directory, strategy.obs_id), file_name)
    fqn = temp_fqn
    if (
        strategy.source_names is not None
        and len(strategy.source_names) > 0
        and os.path.exists(strategy.source_names[0])
        # is there an uncompressed file name?
        and strategy.source_names[0].endswith(file_name)
    ):
        fqn = strategy.source_names[0]
    elif os.path.exists(temp_fqn):
        fqn = temp_fqn
    elif os.path.exists(temp_obs_fqn):
        fqn = temp_obs_fqn
    elif len(strategy.source_names) > 0 and os.path.exists(strategy.source_names[0]):
        # use the compressed file, if it can be found
        fqn = strategy.source_names[0]
    return fqn


class LOTSSPreview:
    """
    Generate a small thumbnail from a previously existing preview image. Override most of the existing
    class, because the original science file is unnecessary.
    """

    def __init__(self, **kwargs):
        self._logger = logging.getLogger(self.__class__.__name__)
        self._release_type = ReleaseType.META
        self._mime_type = 'image/jpeg'
        self._config = kwargs.get('config')
        if self._config is None:
            raise CadcException('Visitor needs a config parameter.')
        self._clients = kwargs.get('clients')
        if self._clients is None or self._clients.data_client is None:
            self._logger.warning('Visitor needs a clients.data_client parameter to store previews.')
        self._strategy = kwargs.get('hierarchy')
        if self._strategy is None:
            raise CadcException('Visitor needs a hierarchy parameter.')
        self._delete_list = []
        # keys are uris, values are lists, where the 0th entry is a file name, and the 1st entry is the artifact type
        self._previews = {}
        self._report = None
        self._hdu_list = None
        self._ext = None
        self._input_fqn = search_for_file(self._strategy, 'preview.jpg', self._strategy.working_directory)
        self._preview_fqn = os.path.join(os.path.dirname(self._input_fqn), self._strategy.prev)
        self._thumb_fqn = os.path.join(os.path.dirname(self._input_fqn), self._strategy.thumb)

    @property
    def report(self):
        return self._report

    def visit(self, observation):
        count = 0
        # preview generation with this algorithm only occurs for the mosaic plane
        if hasattr(self._strategy, 'mosaic_id'):
            product_id = f'{self._strategy.mosaic_id}_mosaic'
            if product_id in observation.planes.keys():
                plane = observation.planes[product_id]
                if not self._strategy.prev_uri in plane.artifacts.keys():
                    self._logger.debug(f'Preview generation for observation {observation.observation_id}, plane {plane.product_id}.')
                    count += self._do_prev(plane, observation.observation_id)
                    self._augment_artifacts(plane)
                    self._delete_list_of_files()
        self._logger.info(f'Changed {count} artifacts during preview augmentation for {observation.observation_id}.')
        self._report = {'artifacts': count}
        return observation

    def generate_plots(self, obs_id):
        count = 0
        if self._strategy._preview_uri:
            self._logger.debug(f'Begin generate_plots for {obs_id} from {self._strategy._preview_uri}')
            http_get(self._strategy._preview_uri, self._input_fqn, self._config.http_get_timeout)
            if os.path.exists(self._input_fqn):
                self._logger.info(f'Retrieved {self._input_fqn}')
                shutil.copy(self._input_fqn, self._preview_fqn)
                self.add_preview(self._strategy.prev_uri, self._strategy.prev, ProductType.PREVIEW)
                count += self._gen_thumbnail()
                if count == 1:
                    self.add_preview(self._strategy.thumb_uri, self._strategy.thumb, ProductType.THUMBNAIL)
        self._logger.debug(f'End generate_plots')
        return count

    def add_preview(self, uri, f_name, product_type, release_type=None, mime_type='image/jpeg'):
        preview_meta = PreviewMeta(f_name, product_type, release_type, mime_type)
        self._previews[uri] = preview_meta

    def add_to_delete(self, fqn):
        self._delete_list.append(fqn)

    def _augment_artifacts(self, plane):
        """Add/update the artifact metadata in the plane."""
        for uri, entry in self._previews.items():
            temp = None
            if uri in plane.artifacts:
                temp = plane.artifacts[uri]
            product_type = entry.product_type
            release_type = self._release_type
            if self._release_type is None:
                release_type = entry.release_type
            if entry.product_type == ProductType.THUMBNAIL:
                fqn = self._thumb_fqn
            else:
                fqn = self._preview_fqn
            plane.artifacts[uri] = get_artifact_metadata(fqn, product_type, release_type, uri, temp)

    def _delete_list_of_files(self):
        """Clean up files on disk after."""
        # cadc_client will be None if executing a ScrapeModify task, so leave the files behind so the user can see
        # them on disk.
        if self._clients is not None and self._clients.data_client is not None:
            for entry in self._delete_list:
                if os.path.exists(entry):
                    self._logger.warning(f'Deleting {entry}')
                    os.unlink(entry)

    def _do_prev(self, plane, obs_id):
        self.generate_plots(obs_id)
        if self._hdu_list is not None:
            # astropy says https://docs.astropy.org/en/stable/io/fits/index.html#working-with-large-files
            self._hdu_list.close()
            del self._hdu_list[self._ext].data
            del self._hdu_list
        self._store_smalls()
        return len(self._previews)

    def _store_smalls(self):
        if self._clients is not None and self._clients.data_client is not None:
            for uri, entry in self._previews.items():
                self._clients.data_client.put(self._strategy.working_directory, uri)

    def _gen_thumbnail(self):
        self._logger.debug(f'Generating thumbnail {self._thumb_fqn}.')
        count = 0
        if os.path.exists(self._preview_fqn):
            # keep import local
            import matplotlib.image as image

            thumb = image.thumbnail(self._preview_fqn, self._thumb_fqn, scale=0.25)
            if thumb is not None:
                count = 1
        else:
            self._logger.warning(f'Could not find {self._preview_fqn} for thumbnail generation.')
        return count

    def _save_figure(self):
        self.add_to_delete(self._preview_fqn)
        count = 1
        self.add_preview(self._strategy.prev_uri, self._strategy.prev, ProductType.PREVIEW, ReleaseType.DATA)
        count += self._gen_thumbnail()
        if count == 2:
            self.add_preview(self._strategy.thumb_uri, self._strategy.thumb, ProductType.THUMBNAIL, ReleaseType.META)
            self.add_to_delete(self._thumb_fqn)
        return count


def visit(observation, **kwargs):
    previewer = LOTSSPreview(**kwargs)
    return previewer.visit(observation)
