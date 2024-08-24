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

"""
This module implements the ObsBlueprint mapping, as well as the workflow entry point that executes the workflow.

"""

import logging

# from os.path import basename, dirname
# from urllib.parse import urlparse
# from astropy import units

from caom2 import CoordCircle2D, DataProductType, PlaneURI, TypedSet, ValueCoord2D
from caom2utils.blueprints import _to_float
from caom2pipe.astro_composable import get_datetime_mjd, get_geocentric_location
from caom2pipe import caom_composable as cc


__all__ = ['mapping_factory']


def _get_float(key, index, metadata):
    value = metadata.get(key)[index]
    return value if isinstance(value, float) else value.item()


def _get_int(key, metadata):
    value = metadata.get(key)
    return value if isinstance(value, int) else value.item()


def _get_int_index(key, index, metadata):
    value = metadata.get(key)[index]
    return value if isinstance(value, int) else value.item()


class DR2MosaicAuxiliaryMapping(cc.TelescopeMapping2):
    def __init__(self, clients, config, dest_uri, hierarchy, observable, observation):
        super().__init__(hierarchy, clients, observable, observation, config)
        self._mosaic_metadata = hierarchy.mosaic_metadata
        self._dest_uri = dest_uri

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model Observation level.

        hard-coded values are from https://science.astron.nl/sdc/astron-data-explorer/data-releases/lotss-dr2/
        """
        self._logger.debug('Begin accumulate_bp.')
        super().accumulate_blueprint(bp)

        # From this page: https://vo.astron.nl/lotss_dr2/q/query_mosaics/info, sidebar item "Created"
        release_date = '2021-08-25T00:00:00.000'
        bp.set('Observation.type', 'OBJECT')
        bp.set('DerivedObservation.members', [])

        bp.set('Observation.algorithm.name', 'mosaic')

        bp.set('Observation.instrument.name', self._mosaic_metadata.get('instid'))

        bp.set('Observation.proposal.id', 'LoTSS')
        bp.set('Observation.proposal.pi', 'T.W. Shimwell')
        bp.set('Observation.proposal.title', 'LOFAR Two-metre Sky Survey')
        bp.set('Observation.proposal.keywords', '_get_observation_proposal_keywords()')

        bp.set('Observation.target.type', 'field')

        telescope_name = 'LOFAR'
        bp.set('Observation.telescope.name', telescope_name)
        x, y, z = get_geocentric_location(telescope_name)
        bp.set('Observation.telescope.geoLocationX', x)
        bp.set('Observation.telescope.geoLocationY', y)
        bp.set('Observation.telescope.geoLocationZ', z)

        bp.set('Plane.metaRelease', release_date)
        bp.set('Plane.dataRelease', release_date)
        bp.set('Plane.calibrationLevel', 4)
        bp.set('Plane.dataProductType', 'image')
        bp.set('Plane.provenance.project', 'LoTSS DR2')
        bp.set('Plane.provenance.producer', 'ASTRON')
        bp.set('Plane.provenance.runID', self._mosaic_metadata.get('data_pid'))
        bp.set('Plane.provenance.lastExecuted', '')
        bp.set('Plane.provenance.reference', self._mosaic_metadata.get('related_products'))

        bp.set('Artifact.productType', self._strategy.get_artifact_product_type(self._dest_uri))
        bp.set('Artifact.releaseType', 'data')
        self._logger.debug('Done accumulate_bp.')

    def update(self):
        """Called to fill multiple CAOM model elements and/or attributes (an n:n relationship between TDM attributes
        and CAOM attributes).
        """
        self._logger.debug('Begin update')
        super().update()
        for plane in self._observation.planes.values():
            # if len(plane.artifacts) >= len(self._storage_name.destination_uris):
            for artifact in plane.artifacts.values():
                for part in artifact.parts.values():
                    for chunk in part.chunks:
                        # no cut-out support for the Temporal Axis
                        chunk.time_axis = None
                        if artifact.uri == self._dest_uri and chunk.naxis == 4 and chunk.polarization is None:
                            # mosaic.fits has only two axes
                            chunk.naxis = 2
                            chunk.energy_axis = None
                        if chunk.naxis == 4 and artifact.uri in [
                            f'astron:LOTSS/{self._strategy.mosaic_id}/mosaic-rms.fits',
                            f'astron:LOTSS/{self._strategy.mosaic_id}/mosaic.pybdsmmask.fits',
                            f'astron:LOTSS/{self._strategy.mosaic_id}/mosaic.resid.fits',
                        ]:
                            # Spectral WCS values from file headers result in negative wavelengths
                            chunk.naxis = 3
        self._logger.debug('End update')
        return self._observation

    def update_time(self):
        pass

    def _get_observation_proposal_keywords(self, ext):
        # values from https://www.aanda.org/articles/aa/full_html/2022/03/aa42484-21/aa42484-21.html
        # default keyword setting splits on whitespace
        temp = set()
        temp.add('surveys')
        temp.add('catalogs')
        temp.add('radio continuum: general')
        temp.add('techniques: image processing')
        return temp

    def _get_provenance_keywords(self, ext):
        # from https://vo.astron.nl/__system__/dc_tables/show/tableinfo/lotss_dr2.mosaics
        d = {
            'C': 'original',
            'F': 'resampled',
            'Z': 'fluxes valid',
            'X': 'not resampled',
            'V': 'for display only',
        }
        temp = set()
        temp.add(d.get(self._mosaic_metadata['pixflags']))
        return temp

    def _get_provenance_name(self, ext):
        temp = self._headers[0].get('ORIGIN')
        result = 'ddf-pipeline'
        if temp:
            result = temp.split()[0]
        return result

    def _get_provenance_version(self, ext):
        temp = self._headers[0].get('ORIGIN')
        result = None
        if temp:
            result = temp.split()[1]
        return result

    def _update_artifact(self, artifact):
        self._logger.debug(f'Begin _update_artifact for {artifact.uri}')
        naxis1 = self._headers[0].get('NAXIS1')
        cdelt1 = self._headers[0].get('CDELT1')
        crval1 = self._headers[0].get('CRVAL1')
        crval2 = self._headers[0].get('CRVAL2')
        if naxis1 and cdelt1 and crval1 and crval2:
            for part in artifact.parts.values():
                for chunk in part.chunks:
                    if chunk.position is not None and chunk.position.axis is not None:
                        # add the circular representation
                        center = ValueCoord2D(
                            coord1=crval1,
                            coord2=crval2,
                        )
                        bounds = CoordCircle2D(center, radius=( naxis1 * abs(cdelt1) / 2.0 ))
                        chunk.position.axis.bounds = bounds
                        self._logger.debug(f'Updated bounds for {self._strategy.file_uri}')

    def _update_plane(self, plane):
        if len(plane.artifacts) > 2:
            mosaic_key = f'{self._strategy.scheme}:{self._strategy.collection}/{self._strategy.mosaic_id}/mosaic.fits'
            if mosaic_key in plane.artifacts.keys():
                mosaic_artifact = plane.artifacts[mosaic_key]
                source_artifact = None
                for artifact in plane.artifacts.values():
                    if artifact.uri != mosaic_key and len(artifact.parts) > 0:
                        source_artifact = artifact
                        break

                if mosaic_artifact and source_artifact:
                    self._logger.info(f'Copy position from {source_artifact.uri} to {mosaic_artifact.uri}')
                    source_chunk = source_artifact.parts['0'].chunks[0]
                    mosaic_chunk = mosaic_artifact.parts['0'].chunks[0]
                    if (
                        source_chunk.position is not None
                        and source_chunk.position.axis is not None
                        and mosaic_chunk.position is not None
                        and mosaic_chunk.position.axis is not None
                    ):
                        mosaic_chunk.position.axis.bounds = source_chunk.position.axis.bounds


class DR2MosaicScience(DR2MosaicAuxiliaryMapping):
    def __init__(self, clients, config, dest_uri, hierarchy, observable, observation):
        super().__init__(clients, config, dest_uri, hierarchy, observable, observation)

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model Observation level."""
        super().accumulate_blueprint(bp)
        bp.set('Plane.provenance.name', '_get_provenance_name()')
        bp.set('Plane.provenance.version', '_get_provenance_version()')
        bp.set('Plane.provenance.keywords', '_get_provenance_keywords()')

        bp.configure_position_axes((1, 2))
        bp.add_attribute('Chunk.position.axis.function.cd11', 'CDELT1')
        bp.set('Chunk.position.axis.function.cd12', 0.0)
        bp.set('Chunk.position.axis.function.cd21', 0.0)
        bp.add_attribute('Chunk.position.axis.function.cd22', 'CDELT2')
        # metadata review comment - ICRS
        bp.set('Chunk.position.coordsys', 'ICRS')
        # 6"
        bp.set('Chunk.position.resolution', 0.001666)

        self._logger.debug('Done accumulate_bp.')


class DR2MosaicScienceLow(DR2MosaicScience):
    def __init__(self, clients, config, dest_uri, hierarchy, observable, observation):
        super().__init__(clients, config, dest_uri, hierarchy, observable, observation)

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model Observation level."""
        super().accumulate_blueprint(bp)
        # 20"
        bp.set('Chunk.position.resolution', 0.00555)
        bp.configure_energy_axis(4)
        # TODO how to translate "2 channels per 0.195 MHz subband" to resolution?
        bp.set('Chunk.energy.specsys', 'TOPOCENT')
        bp.set('Chunk.energy.bandpassName', self._mosaic_metadata['bandpassid'])
        bp.set('Chunk.energy.axis.axis.ctype', 'WAVE')
        bp.set('Chunk.energy.axis.axis.cunit', self._mosaic_metadata['bandpassunit'])
        bp.set('Chunk.energy.axis.range.start.pix', 0.5)
        bp.set('Chunk.energy.axis.range.start.val', self._mosaic_metadata['bandpasshi'])
        bp.set('Chunk.energy.axis.range.end.pix', 1.5)
        bp.set('Chunk.energy.axis.range.end.val', self._mosaic_metadata['bandpasslo'])
        self._logger.debug('Done accumulate_bp.')


class DR2Mosaic(DR2MosaicAuxiliaryMapping):
    def __init__(self, clients, config, dest_uri, hierarchy, observable, observation):
        super().__init__(clients, config, dest_uri, hierarchy, observable, observation)

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model Observation level."""
        super().accumulate_blueprint(bp)

        bp.set('Artifact.contentType', 'application/fits')
        bp.set('Artifact.contentLength', _get_int('accsize', self._mosaic_metadata))

        bp.configure_position_axes((1, 2))
        bp.set('Chunk.position.coordsys', self._mosaic_metadata.get('refframe'))
        bp.set('Chunk.position.equinox', self._mosaic_metadata.get('wcs_equinox'))
        # bp.set('Chunk.position.axis.axis1.ctype', self._mosaic_metadata['wcs_projection'])
        bp.set('Chunk.position.axis.axis1.ctype', 'RA---SIN')
        bp.set('Chunk.position.axis.axis1.cunit', 'deg')
        # bp.set('Chunk.position.axis.axis2.ctype', self._mosaic_metadata['wcs_projection'])
        bp.set('Chunk.position.axis.axis2.ctype', 'DEC--SIN')
        bp.set('Chunk.position.axis.axis2.cunit', 'deg')
        bp.set('Chunk.position.axis.function.cd11', _get_float('wcs_cdmatrix', 0, self._mosaic_metadata))
        bp.set('Chunk.position.axis.function.cd12', _get_float('wcs_cdmatrix', 1, self._mosaic_metadata))
        bp.set('Chunk.position.axis.function.cd21', _get_float('wcs_cdmatrix', 2, self._mosaic_metadata))
        bp.set('Chunk.position.axis.function.cd22', _get_float('wcs_cdmatrix', 3, self._mosaic_metadata))
        bp.set('Chunk.position.axis.function.dimension.naxis1', _get_int_index('pixelsize', 0, self._mosaic_metadata))
        bp.set('Chunk.position.axis.function.dimension.naxis2', _get_int_index('pixelsize', 1, self._mosaic_metadata))
        bp.set('Chunk.position.axis.function.refCoord.coord1.pix',  _get_float('wcs_refpixel', 0, self._mosaic_metadata))
        bp.set('Chunk.position.axis.function.refCoord.coord1.val',  _get_float('wcs_refvalues', 0, self._mosaic_metadata))
        bp.set('Chunk.position.axis.function.refCoord.coord2.pix',  _get_float('wcs_refpixel', 1, self._mosaic_metadata))
        bp.set('Chunk.position.axis.function.refCoord.coord2.val',  _get_float('wcs_refvalues', 1, self._mosaic_metadata))
        # 6"
        bp.set('Chunk.position.resolution', 0.001666)

        bp.configure_energy_axis(4)
        # TODO how to translate "2 channels per 0.195 MHz subband" to resolution?
        bp.set('Chunk.energy.specsys', 'TOPOCENT')
        bp.set('Chunk.energy.bandpassName', self._mosaic_metadata['bandpassid'])
        bp.set('Chunk.energy.axis.axis.ctype', 'WAVE')
        bp.set('Chunk.energy.axis.axis.cunit', self._mosaic_metadata['bandpassunit'])
        bp.set('Chunk.energy.axis.range.start.pix', 0.5)
        bp.set('Chunk.energy.axis.range.start.val', self._mosaic_metadata['bandpasshi'])
        bp.set('Chunk.energy.axis.range.end.pix', 1.5)
        bp.set('Chunk.energy.axis.range.end.val', self._mosaic_metadata['bandpasslo'])

        self._logger.debug('Done accumulate_bp.')

    def _update_plane(self, plane):
        self._logger.debug('Begin _update_plane')
        if len(plane.artifacts) > 1:
            mosaic_artifact = None
            source_artifact = None
            for artifact in plane.artifacts.values():
                if artifact.uri == self._strategy.file_name:
                    mosaic_artifact = artifact
                else:
                    if len(artifact.parts) > 0:
                        source_artifact = artifact
                if mosaic_artifact and source_artifact:
                    break

            if mosaic_artifact and source_artifact:
                self._logger.info(f'Copy position from {source_artifact.uri} to {mosaic_artifact.uri}')
                source_chunk = source_artifact.parts['0'].chunks[0]
                mosaic_chunk = mosaic_artifact.parts['0'].chunks[0]
                if (
                    source_chunk.position is not None
                    and source_chunk.position.axis is not None
                    and mosaic_chunk.position is not None
                    and mosaic_chunk.position.axis is not None
                ):
                    mosaic_chunk.position.axis.bounds = source_chunk.position.axis.bounds
        self._logger.debug('End _update_plane')


class DR2MosaicSciencePolarization(DR2MosaicScience):
    def __init__(self, clients, config, dest_uri, hierarchy, observable, observation):
        super().__init__(clients, config, dest_uri, hierarchy, observable, observation)

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model Observation level."""
        super().accumulate_blueprint(bp)
        bp.set('Plane.provenance.name', '_get_provenance_name()')
        bp.set('Plane.provenance.version', '_get_provenance_version()')
        bp.set('Plane.provenance.keywords', '_get_provenance_keywords()')

        # if I use the 4th axis values from the files mosaic.resid.fits, mosaic-rms.fits, and mosaic.pybdsmmask.fits
        # the Frequency ends up negative at the Plane level.
        #
        # Leave those values out for now
        #
        # bp.configure_energy_axis(4)
        # # # TODO how to translate "2 channels per 0.195 MHz subband" to resolution?
        # bp.set('Chunk.energy.bandpassName', self._mosaic_metadata['bandpassid'])

        bp.configure_polarization_axis(3)
        self._logger.debug('Done accumulate_bp.')


class DR2MosaicSciencePolarizationLow(DR2MosaicSciencePolarization):
    def __init__(self, clients, config, dest_uri, hierarchy, observable, observation):
        super().__init__(clients, config, dest_uri, hierarchy, observable, observation)

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model Observation level."""
        super().accumulate_blueprint(bp)
        # 20"
        bp.set('Chunk.position.resolution', 0.00555)
        self._logger.debug('Done accumulate_bp.')


class DR2Raw(cc.TelescopeMapping2):
    def __init__(self, clients, config, dest_uri, hierarchy, observable, observation):
        super().__init__(hierarchy, clients, observable, observation, config)
        self._dest_uri = dest_uri

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model Observation level.

        hard-coded values are from https://science.astron.nl/sdc/astron-data-explorer/data-releases/lotss-dr2/
        """
        self._logger.debug('Begin accumulate_bp.')
        super().accumulate_blueprint(bp)

        release_date = self._strategy.metadata.release_date
        bp.set('Observation.type', 'OBJECT')
        bp.set('Observation.target.type', 'field')

        telescope_name = 'LOFAR'
        bp.set('Observation.telescope.name', telescope_name)
        x, y, z = get_geocentric_location(telescope_name)
        bp.set('Observation.telescope.geoLocationX', x)
        bp.set('Observation.telescope.geoLocationY', y)
        bp.set('Observation.telescope.geoLocationZ', z)

        bp.set('Plane.metaRelease', release_date)
        bp.set('Plane.dataRelease', release_date)
        bp.set('Plane.calibrationLevel', 1)
        bp.set('Plane.dataProductType', 'measurements')
        bp.set('Plane.provenance.name', self._strategy.metadata.project_name)
        bp.set('Plane.provenance.producer', self._strategy.metadata.creator)
        bp.set('Plane.provenance.lastExecuted', '')
        # TODO
        # bp.set('Plane.provenance.reference', self._mosaic_metadata.get('related_products'))

        bp.set('Artifact.productType', self._strategy.get_artifact_product_type(self._dest_uri))
        bp.set('Artifact.releaseType', 'data')
        bp.set('Artifact.contentChecksum', self._strategy.metadata.content_checksum)
        bp.set('Artifact.contentType', self._strategy.metadata.content_type)

        bp.configure_position_axes((1, 2))
        bp.set('Chunk.position.axis.axis1.ctype', 'RA---SIN')
        bp.set('Chunk.position.axis.axis1.cunit', 'deg')
        bp.set('Chunk.position.axis.axis2.ctype', 'DEC--SIN')
        bp.set('Chunk.position.axis.axis2.cunit', 'deg')
        bp.add_attribute('Chunk.position.axis.function.cd11', 'CDELT1')
        bp.set('Chunk.position.axis.function.cd12', 0.0)
        bp.set('Chunk.position.axis.function.cd21', 0.0)
        bp.add_attribute('Chunk.position.axis.function.cd22', 'CDELT2')
        bp.set('Chunk.position.coordsys', 'ICRS')

        bp.configure_energy_axis(3)
        bp.set('Chunk.energy.axis.axis.ctype', 'WAVE')
        bp.set('Chunk.energy.axis.axis.cunit', 'MHz')
        # Central Frequency
        # Channel Width
        # self._logger.error(self._strategy.metadata)
        central_frequency = _to_float(self._strategy.metadata.central_frequency)  # MHz
        channel_width = _to_float(self._strategy.metadata.channel_width) / 10e6  # Central Frequency in Hz
        bp.set('Chunk.energy.axis.range.start.pix', 0.5)
        bp.set('Chunk.energy.axis.range.start.val', central_frequency - (channel_width / 2.0))
        bp.set('Chunk.energy.axis.range.end.pix', 1.5)
        bp.set('Chunk.energy.axis.range.end.val',  central_frequency + (channel_width / 2.0))
        bp.set('Chunk.energy.resolvingPower', self._strategy.metadata.integration_interval)
        # TODO - guess
        bp.set('Chunk.energy.specsys', 'TOPOCENT')

        # 22-08-24 Aida Ahmadi
        # There is a caveat that in interferometry the integration time is defined as something else (usually on the
        # order of ~1 second → 'Integration Interval' on the LTA). This is a limitation of having to adapt terms
        # from optical astronomy that have been used as standard to radio and interferometry needs. On the CAOM
        # side, the values should also be updated in the Chunks. From looking at the Chunk.time explanation here,
        # this should be ingested under 'exposure' rather than 'resolution' and 'cdelt' where it is currently
        # ingested.
        bp.configure_time_axis(4)
        bp.set('Chunk.time.axis.axis.ctype', 'TIME')
        bp.set('Chunk.time.axis.axis.cunit', 'd')
        bp.set('Chunk.time.axis.range.start.pix', 0.5)
        start_mjd = get_datetime_mjd(self._strategy.metadata.start_time)
        if start_mjd:
            bp.set('Chunk.time.axis.range.start.val', start_mjd.value)
        bp.set('Chunk.time.axis.range.end.pix', 1.5)
        end_mjd = get_datetime_mjd(self._strategy.metadata.end_time)
        if end_mjd:
            bp.set('Chunk.time.axis.range.end.val', end_mjd.value)
        # exposure units are 's'
        bp.set('Chunk.time.exposure', _to_float(self._strategy.metadata.duration))
        bp.set('Chunk.time.timesys', 'UTC')

        self._logger.debug('Done accumulate_bp.')

    def update(self):
        """Called to fill multiple CAOM model elements and/or attributes (an n:n relationship between TDM attributes
        and CAOM attributes).
        """
        super().update()
        provenance_plane_uris = TypedSet(PlaneURI,)
        for plane in self._observation.planes.values():
            if plane.data_product_type == DataProductType.MEASUREMENTS:
                _, plane_uri = cc.make_plane_uri(
                    self._strategy.obs_id, self._strategy.product_id, self._strategy.collection
                )
                provenance_plane_uris.add(plane_uri)

        if len(provenance_plane_uris) > 0:
            for plane in self._observation.planes.values():
                if plane.data_product_type != DataProductType.MEASUREMENTS:
                    plane.provenance.inputs.update(provenance_plane_uris)
                for artifact in plane.artifacts.values():
                    for part in artifact.parts.values():
                        for chunk in part.chunks:
                            # no cut-out support for any axes
                            # chunk.position_axis_1 = 1
                            # chunk.position_axis_2 = 2
                            chunk.time_axis = None
                            chunk.energy_axis = None
        return self._observation

    def _get_observation_proposal_keywords(self, ext):
        # values from https://www.aanda.org/articles/aa/full_html/2022/03/aa42484-21/aa42484-21.html
        # default keyword setting splits on whitespace
        temp = set()
        temp.add('surveys')
        temp.add('catalogs')
        temp.add('radio continuum: general')
        temp.add('techniques: image processing')
        return temp

    def _update_artifact(self, artifact):
        self._logger.debug(f'Begin _update_artifact for {artifact.uri}')
        # naxis1 = self._strategy.metadata.get('TODO')
        # current values are copied from the processed records
        naxis1 = 2978
        # cdelt1 = self._strategy.metadata.get('TODO')
        cdelt1 = -0.00041666666666666
        crval1 = self._strategy.metadata.ra
        crval2 = self._strategy.metadata.dec
        if naxis1 and cdelt1 and crval1 and crval2:
            for part in artifact.parts.values():
                for chunk in part.chunks:
                    if chunk.position is not None and chunk.position.axis is not None:
                        # add the circular representation
                        center = ValueCoord2D(coord1=crval1, coord2=crval2)
                        bounds = CoordCircle2D(center, radius=(naxis1 * abs(cdelt1) / 2.0 ))
                        chunk.position.axis.bounds = bounds
                        self._logger.debug(f'Updated bounds for {self._strategy.file_uri}')
                    # TODO like this until the Spatial WCS handling is organized
                    chunk.naxis = None
                    chunk.position_axis_1 = None
                    chunk.position_axis_2 = None

    def update_time(self):
        """
        22-08-24 Aida Ahmadi

        The values should come from the awe query for the measurement set(s) associated with the FITS image but on
        the 'Interferometric Data' level and not the 'Averaging Pipeline' level (i.e. here:
        https://lta.lofar.eu/Lofar?project=ALL&mode=query_result_page&product=CorrelatedDataProduct&pipeline_object_id=740B81E0588B278FE053144A17ACFF94
        and not
        https://lta.lofar.eu/Lofar?project=ALL&product=all_observation_pipeline&mode=query_result_page_user&ObservationId=664568).

        The latter (Averaging Pipeline) start/end times refer to how long the averaging pipeline took to run.

        Use the start/end times of the measurement sets associated with mosaic and mosaic_low images and to get
        the durations also from the awe query.

        In some cases, to reach the needed sensitivity (usually for low declination sources) a different observing
        strategy was adopted where a few blocks of shorter duration observations were made. This means that to produce
        the FITS images for low declination sources, more than one group of measurement sets have been used. An
        example is target P005+21, where at the bottom of the list on the ASTRON VO page
        (https://vo.astron.nl/lotss_dr2/q/dlmosaic/dlmeta?ID=ivo%3A//astron.nl/~%3FLoTSS-DR2/P005%2B21) you can see
        two progenitor links. One that goes to SAS ID 727372
        (https://lta.lofar.eu/Lofar?project=ALL&product=all_observation_pipeline&mode=query_result_page_user&ObservationId=727372)
        and one that goes to SAS ID 734329
        (https://lta.lofar.eu/Lofar?project=ALL&product=all_observation_pipeline&mode=query_result_page_user&ObservationId=734329).

        Now in such a case, observations can even be years apart. Use the earliest start date and the latest end date from the
        LTA awe query at the 'Interferometric Data' level (i.e.
        https://lta.lofar.eu/Lofar?project=ALL&mode=query_result_page&product=CorrelatedDataProduct&pipeline_object_id=90C703E9A9705084E053144A17AC6B3C
        and
        https://lta.lofar.eu/Lofar?project=ALL&mode=query_result_page&product=CorrelatedDataProduct&pipeline_object_id=8CEF25CF085D6A53E053164A17AC7525)

        AND to include in one of the CADC UI columns and in the CAOM model the total duration which should be the sum
        of the two durations reported for each of the processes: 28799.0+14399.0 (in seconds, so ~12 hours).
        """
        if len(self._observation.planes) >= 3:
            from caom2 import Axis, CoordAxis1D, CoordBounds1D, CoordRange1D, DataProductType, RefCoord, TemporalWCS, TypedList

            # assume that the MS Artifacts all have the same TemporalWCS, so the first instance is sufficient to build
            # the collection
            copy_from_chunks = []
            for plane in self._observation.planes.values():
                if plane.data_product_type == DataProductType.MEASUREMENTS:
                    for artifact in plane.artifacts.values():
                        for part in artifact.parts.values():
                            for chunk in part.chunks:
                                if chunk.time and chunk.time.axis and chunk.time.axis.range:
                                    copy_from_chunks.append(chunk)
                                    break
                            break
                        break

            if len(copy_from_chunks) > 0:
                samples = TypedList(CoordRange1D,)
                exposure = 0.0
                for chunk in copy_from_chunks:
                    exposure += chunk.time.exposure
                    start_ref_coord = RefCoord(pix=0.5, val=chunk.time.axis.range.start.val)
                    end_ref_coord = RefCoord(pix=1.5, val=chunk.time.axis.range.end.val)
                    sample = CoordRange1D(start_ref_coord, end_ref_coord)
                    samples.append(sample)

                bounds = CoordBounds1D(samples=samples)
                axis = CoordAxis1D(Axis(ctype='TIME', cunit='d'), bounds=bounds)
                t = TemporalWCS(axis=axis, timesys='UTC', exposure=exposure)

                for plane in self._observation.planes.values():
                    if plane.data_product_type != DataProductType.MEASUREMENTS:
                        for artifact in plane.artifacts.values():
                            for part in artifact.parts.values():
                                for chunk in part.chunks:
                                    chunk.time = t


def mapping_factory(clients, config, dest_uri, hierarchy, observable, observation):
    result = None
    if dest_uri.endswith('MS') or dest_uri.endswith('tar'):
        result = DR2Raw(clients, config, dest_uri, hierarchy, observable, observation)
    elif dest_uri == f'{hierarchy.scheme}:{hierarchy.collection}/{hierarchy.mosaic_id}/mosaic-blanked.fits':
        result = DR2Mosaic(clients, config, dest_uri, hierarchy, observable, observation)
    else:
        if hierarchy.product_id.endswith('low'):
            if hierarchy.metadata[0].get('NAXIS') == 2:
                result = DR2MosaicScienceLow(clients, config, dest_uri, hierarchy, observable, observation)
            else:
                result = DR2MosaicSciencePolarizationLow(clients, config, dest_uri, hierarchy, observable, observation)
        else:
            if hierarchy.metadata[0].get('NAXIS') == 2:
                result = DR2MosaicScience(clients, config, dest_uri, hierarchy, observable, observation)
            else:
                result = DR2MosaicSciencePolarization(clients, config, dest_uri, hierarchy, observable, observation)
    logging.debug(f'Created {result.__class__.__name__} for {dest_uri}')
    return result
