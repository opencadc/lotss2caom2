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

import logging

from caom2 import Observation, DataProductType
from caom2pipe import caom_composable as cc
from caom2pipe import manage_composable as mc


__all__ = ['visit']


def visit(observation, **kwargs):
    assert observation is not None, 'Input parameter must have a value.'
    assert isinstance(observation, Observation), 'Input parameter must be an Observation'

    working_dir = kwargs.get('working_directory', './')
    hierarchy = kwargs.get('hierarchy')
    if hierarchy is None:
        raise mc.CadcException(f'No hierarchy provided to visitor for obs {observation.observation_id}.')
    log_file_directory = kwargs.get('log_file_directory')

    logging.info(f'Begin footprint finding for {hierarchy.get_file_fqn(working_dir)}.')
    count = 0
    original_chunk = None
    for plane in observation.planes.values():
        if plane.data_product_type != DataProductType.MEASUREMENTS:
            for artifact in plane.artifacts.values():
                if artifact.uri.endswith('mosaic-blanked.fits'):
                    for part in artifact.parts.values():
                        for chunk in part.chunks:
                            # -t 10 provides a margin of up to 10 pixels
                            cc.exec_footprintfinder(
                                chunk,
                                hierarchy.get_file_fqn(working_dir),
                                log_file_directory,
                                hierarchy.file_id,
                                '-t 10',
                            )
                            original_chunk = chunk
                            count += 1
                            break
                        if count == 1:
                            break  # part
                if count == 1:
                    break  # artifact
        if count == 1:
            break  # plane

    if original_chunk:
        for plane in observation.planes.values():
            if plane.data_product_type != DataProductType.MEASUREMENTS:
                for artifact in plane.artifacts.values():
                    if not artifact.uri.endswith('/mosaic-blanked.fits'):
                        for part in artifact.parts.values():
                            for chunk in part.chunks:
                                chunk.position = original_chunk.position
                                count += 1

    logging.info(f'Completed footprint augmentation. Changed {count} artifacts.')
    return observation
