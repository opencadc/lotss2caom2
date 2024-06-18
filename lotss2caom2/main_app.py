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
This module implements the ObsBlueprint mapping, as well as the workflow
entry point that executes the workflow.
"""

import logging

# from os.path import basename, dirname
# from urllib.parse import urlparse
# from astropy import units

from caom2 import CoordCircle2D, ProductType, ValueCoord2D
from caom2utils.blueprints import _to_float
from caom2pipe.astro_composable import get_datetime_mjd, get_geocentric_location
from caom2pipe import caom_composable as cc


__all__ = ['mapping_factory']


class DR2MosaicAuxiliaryMapping(cc.TelescopeMapping2):
    def __init__(self, strategy, clients, observable, observation, config, dest_uri):
        super().__init__(strategy, clients, observable, observation, config)
        self._mosaic_metadata = strategy.mosaic_metadata
        self._dest_uri = dest_uri

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model Observation level.

        hard-coded values are from https://science.astron.nl/sdc/astron-data-explorer/data-releases/lotss-dr2/
        """
        self._logger.debug('Begin accumulate_bp.')
        super().accumulate_blueprint(bp)

        release_date = '2023-01-01T00:00:00.000'
        bp.set('Observation.metaRelease', release_date)
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

        bp.configure_time_axis(5)
        bp.set('Chunk.time.axis.axis.ctype', 'TIME')
        bp.set('Chunk.time.axis.axis.cunit', 'd')
        bp.set('Chunk.time.axis.function.naxis', 1)
        # TODO - query to determine whether the observation was 8 or 16 hours
        bp.set('Chunk.time.axis.function.delta', 8 / 24.0)
        bp.set('Chunk.time.axis.function.refCoord.pix', 0.5)
        bp.set('Chunk.time.axis.function.refCoord.val', self._mosaic_metadata.get('dateobs'))
        # resolution units are 'd'
        bp.set('Chunk.time.resolution', 8 / 24.0)
        bp.set('Chunk.time.timesys', 'UTC')

        self._logger.debug('Done accumulate_bp.')

    def update(self):
        """Called to fill multiple CAOM model elements and/or attributes (an n:n relationship between TDM attributes
        and CAOM attributes).
        """
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
                    self._logger.error(f'Copy position from {source_artifact.uri} to {mosaic_artifact.uri}')
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
    def __init__(self, strategy, clients, observable, observation, config, dest_uri):
        super().__init__(strategy, clients, observable, observation, config, dest_uri)

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
    def __init__(self, strategy, clients, observable, observation, config, dest_uri):
        super().__init__(strategy, clients, observable, observation, config, dest_uri)

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
    def __init__(self, strategy, clients, observable, observation, config, dest_uri):
        super().__init__(strategy, clients, observable, observation, config, dest_uri)

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model Observation level."""
        super().accumulate_blueprint(bp)

        bp.set('Artifact.contentType', 'application/fits')
        bp.set('Artifact.contentLength', self._mosaic_metadata.get('accsize').item())

        bp.configure_position_axes((1, 2))
        bp.set('Chunk.position.coordsys', self._mosaic_metadata.get('refframe'))
        bp.set('Chunk.position.equinox', self._mosaic_metadata.get('wcs_equinox'))
        # bp.set('Chunk.position.axis.axis1.ctype', self._mosaic_metadata['wcs_projection'])
        bp.set('Chunk.position.axis.axis1.ctype', 'RA---SIN')
        bp.set('Chunk.position.axis.axis1.cunit', 'deg')
        # bp.set('Chunk.position.axis.axis2.ctype', self._mosaic_metadata['wcs_projection'])
        bp.set('Chunk.position.axis.axis2.ctype', 'DEC--SIN')
        bp.set('Chunk.position.axis.axis2.cunit', 'deg')
        cd_matrix = self._mosaic_metadata.get('wcs_cdmatrix')
        bp.set('Chunk.position.axis.function.cd11', cd_matrix[0].item())
        bp.set('Chunk.position.axis.function.cd12', cd_matrix[1].item())
        bp.set('Chunk.position.axis.function.cd21', cd_matrix[2].item())
        bp.set('Chunk.position.axis.function.cd22', cd_matrix[3].item())
        pixel_size = self._mosaic_metadata.get('pixelsize')
        bp.set('Chunk.position.axis.function.dimension.naxis1', pixel_size[0].item())
        bp.set('Chunk.position.axis.function.dimension.naxis2', pixel_size[1].item())
        ref_pixel = self._mosaic_metadata.get('wcs_refpixel')
        ref_values = self._mosaic_metadata.get('wcs_refvalues')
        bp.set('Chunk.position.axis.function.refCoord.coord1.pix', ref_pixel[0].item())
        bp.set('Chunk.position.axis.function.refCoord.coord1.val', ref_values[0].item())
        bp.set('Chunk.position.axis.function.refCoord.coord2.pix', ref_pixel[1].item())
        bp.set('Chunk.position.axis.function.refCoord.coord2.val', ref_values[1].item())
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
                self._logger.error(f'Copy position from {source_artifact.uri} to {mosaic_artifact.uri}')
                source_chunk = source_artifact.parts['0'].chunks[0]
                mosaic_chunk = mosaic_artifact.parts['0'].chunks[0]
                if (
                    source_chunk.position is not None 
                    and source_chunk.position.axis is not None 
                    and mosaic_chunk.position is not None 
                    and mosaic_chunk.position.axis is not None 
                ):
                    mosaic_chunk.position.axis.bounds = source_chunk.position.axis.bounds


class DR2MosaicSciencePolarization(DR2MosaicScience):
    def __init__(self, strategy, clients, observable, observation, config, dest_uri):
        super().__init__(strategy, clients, observable, observation, config, dest_uri)

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
    def __init__(self, strategy, clients, observable, observation, config, dest_uri):
        super().__init__(strategy, clients, observable, observation, config, dest_uri)

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model Observation level."""
        super().accumulate_blueprint(bp)
        # 20"
        bp.set('Chunk.position.resolution', 0.00555)
        self._logger.debug('Done accumulate_bp.')


class DR2Raw(cc.TelescopeMapping2):
    def __init__(self, strategy, clients, observable, observation, config, dest_uri):
        super().__init__(strategy, clients, observable, observation, config)
        self._dest_uri = dest_uri

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model Observation level.

        hard-coded values are from https://science.astron.nl/sdc/astron-data-explorer/data-releases/lotss-dr2/
        """
        self._logger.debug('Begin accumulate_bp.')
        super().accumulate_blueprint(bp)

        release_date = '2023-01-01T00:00:00.000'
        bp.set('Observation.metaRelease', release_date)
        bp.set('Observation.type', 'OBJECT')
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
        bp.set('Plane.calibrationLevel', 1)
        bp.set('Plane.dataProductType', 'measurements')
        bp.set('Plane.provenance.project', 'LoTSS DR2')
        bp.set('Plane.provenance.producer', self._strategy.metadata.get('Creator'))
        bp.set('Plane.provenance.lastExecuted', '')
        # TODO
        # bp.set('Plane.provenance.reference', self._mosaic_metadata.get('related_products'))

        bp.set('Artifact.productType', self._strategy.get_artifact_product_type(self._dest_uri))
        bp.set('Artifact.releaseType', 'data')

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
        central_frequency = _to_float(self._strategy.metadata.get('Central Frequency [MHz]'))
        channel_width = _to_float(self._strategy.metadata.get('Channel Width [Hz]')) / 1*10e6  # Central Frequency units
        bp.set('Chunk.energy.axis.range.start.pix', 0.5)
        bp.set('Chunk.energy.axis.range.start.val', central_frequency - (channel_width / 2.0))
        bp.set('Chunk.energy.axis.range.end.pix', 1.5)
        bp.set('Chunk.energy.axis.range.end.val',  central_frequency + (channel_width / 2.0))
        bp.set('Chunk.energy.resolvingPower', self._strategy.metadata.get('Integration Interval [s]'))
        # TODO - guess
        bp.set('Chunk.energy.specsys', 'TOPOCENT')

        bp.configure_time_axis(4)
        bp.set('Chunk.time.axis.axis.ctype', 'TIME')
        bp.set('Chunk.time.axis.axis.cunit', 'd')
        bp.set('Chunk.time.axis.range.start.pix', 0.5)
        start_mjd = get_datetime_mjd(self._strategy.metadata.get('Start Time'))
        if start_mjd:
            bp.set('Chunk.time.axis.range.start.val', start_mjd.value)
        bp.set('Chunk.time.axis.range.end.pix', 1.5)
        end_mjd = get_datetime_mjd(self._strategy.metadata.get('End Time'))
        if end_mjd:
            bp.set('Chunk.time.axis.range.end.val', end_mjd.value)
        bp.set('Chunk.time.resolution', _to_float(self._strategy.metadata.get('Duration [s]')) / (24.0 * 3600.0) )
        bp.set('Chunk.time.timesys', 'UTC')

        self._logger.debug('Done accumulate_bp.')

    def update(self):
        """Called to fill multiple CAOM model elements and/or attributes (an n:n relationship between TDM attributes
        and CAOM attributes).
        """
        super().update()
        for plane in self._observation.planes.values():
            for artifact in plane.artifacts.values():
                for part in artifact.parts.values():
                    for chunk in part.chunks:
                        # no cut-out support for any axes
                        chunk.position_axis_1 = 1
                        chunk.position_axis_2 = 2
                        chunk.time_axis = None
                        chunk.energy_axis = None
                        chunk.naxis = 0
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
        naxis1 = self._strategy.metadata.get('TODO')
        cdelt1 = self._strategy.metadata.get('TODO')
        crval1 = self._strategy.metadata.get('Right Ascension [degrees]')
        crval2 = self._strategy.metadata.get('Declination [degrees]')
        if naxis1 and cdelt1 and crval1 and crval2:
            for part in artifact.parts.values():
                for chunk in part.chunks:
                    if chunk.position is not None and chunk.position.axis is not None:
                        # add the circular representation
                        center = ValueCoord2D(coord1=crval1, coord2=crval2)
                        bounds = CoordCircle2D(center, radius=(naxis1 * abs(cdelt1) / 2.0 ))
                        chunk.position.axis.bounds = bounds
                        self._logger.debug(f'Updated bounds for {self._strategy.file_uri}')


def mapping_factory(strategy, clients, observable, observation, config, dest_uri):
    # logging.error('\n'.join(ii for ii in storage_name._headers.keys()))
    # logging.error(dest_uri)
    result = None
    if dest_uri.endswith('MS'):
        result = DR2Raw(strategy, clients, observable, observation, config, dest_uri)
    elif dest_uri == f'{strategy.scheme}:{strategy.collection}/{strategy.mosaic_id}/mosaic.fits':
        result = DR2Mosaic(strategy, clients, observable, observation, config, dest_uri)
    else:
        if strategy.product_id.endswith('low'):
            if strategy.metadata[0].get('NAXIS') == 2:
                result = DR2MosaicScienceLow(strategy, clients, observable, observation, config, dest_uri)
            else:
                result = DR2MosaicSciencePolarizationLow(strategy, clients, observable, observation, config, dest_uri)
        else:
            if strategy.metadata[0].get('NAXIS') == 2:
                result = DR2MosaicScience(strategy, clients, observable, observation, config, dest_uri)
            else:
                result = DR2MosaicSciencePolarization(strategy, clients, observable, observation, config, dest_uri)
    logging.debug(f'Created {result.__class__.__name__} for {dest_uri}')
    return result
