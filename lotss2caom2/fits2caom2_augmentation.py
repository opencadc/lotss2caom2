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
#  : 4 $
#
# ***********************************************************************
#

import os
import traceback

from caom2 import ProductType, Algorithm, SimpleObservation, DerivedObservation
from caom2utils import ContentParser, FitsParser
from caom2utils.parsers import Caom2Exception
from caom2pipe import caom_composable as cc
from lotss2caom2 import main_app


__all__ = ['LoTSSFits2caom2Visitor']


class LoTSSFits2caom2Visitor(cc.Fits2caom2Visitor):
    def __init__(self, observation, **kwargs):
        super().__init__(observation, **kwargs)
        self._strategy = kwargs.get('hierarchy')

    def _get_mapping(self, dest_uri):
        return main_app.mapping_factory(
            self._strategy,
            self._clients,
            self._observable,
            self._observation,
            self._config,
            dest_uri,
        )

    def _get_parser(self, blueprint, uri):
        if (
            hasattr(uri, 'mosaic')
            and uri == f'{self._strategy.scheme}:{self._strategy.collection}/{self._strategy.mosaic_id}/mosaic.fits'
        ):
            parser = ContentParser(blueprint, uri)
        else:
            parser = FitsParser(self._strategy.metadata, blueprint, uri)
        self._logger.debug(f'Creating {parser.__class__.__name__} for {uri}')
        return parser

    def visit(self):
        self._logger.debug('Begin visit')
        try:
            for uri in self._strategy.destination_uris:
                self._logger.info(f'Build observation for {uri}')
                telescope_data = self._get_mapping(uri)
                if telescope_data is None:
                    self._logger.info(f'Ignoring {uri} because there is no TelescopeMapping.')
                    continue
                blueprint = self._get_blueprint(telescope_data)
                telescope_data.accumulate_blueprint(blueprint)
                if self._config.dump_blueprint and self._config.log_to_file:
                    with open(f'{self._config.log_file_directory}/{os.path.basename(uri)}.bp', 'w') as f:
                        f.write(blueprint.__str__())
                parser = self._get_parser(blueprint, uri)

                if self._observation is None:
                    if blueprint._get('DerivedObservation.members') is None:
                        self._logger.debug('Build a SimpleObservation')
                        self._observation = SimpleObservation(
                            collection=self._strategy.collection,
                            observation_id=self._strategy.obs_id,
                            algorithm=Algorithm('exposure'),
                        )
                    else:
                        self._logger.debug('Build a DerivedObservation')
                        algorithm_name =(
                            'composite'
                            if blueprint._get('Observation.algorithm.name') == 'exposure'
                            else blueprint._get('Observation.algorithm.name')
                        )
                        self._observation = DerivedObservation(
                            collection=self._strategy.collection,
                            observation_id=self._strategy.obs_id,
                            algorithm=Algorithm(algorithm_name),
                        )
                    telescope_data.observation = self._observation
                parser.augment_observation(
                    observation=self._observation,
                    artifact_uri=uri,
                    product_id=self._strategy.product_id,
                )

                # file_info = self._metadata_reader.file_info.get(uri)
                self._observation = telescope_data.update()
        except Caom2Exception as e:
            self._logger.debug(traceback.format_exc())
            self._logger.warning(
                f'CAOM2 record creation failed for {self._strategy.obs_id}:{self._strategy.file_name} with {e}'
            )
            self._observation = None

        self._logger.debug('End visit')
        return self._observation


def visit(observation, **kwargs):
    return LoTSSFits2caom2Visitor(observation, **kwargs).visit()
