import sys
from astropy import units
from caom2_24 import ObservationReader as or_24
from caom2_24 import Circle as Circle_24
from caom2_24 import Position as Position_24
from caom2 import ObservationWriter as ow_25
from caom2 import Algorithm, Artifact, CalibrationLevel, Circle, DataLinkSemantics, DataProductType
from caom2 import Dimension2D, Energy, EnergyBand, Instrument
from caom2 import MultiShape, ObservationIntentType, Plane, Point, Polarization, PolarizationState
from caom2 import Polygon, Position
from caom2 import Proposal, Provenance
from caom2 import ReleaseType, SimpleObservation, Target, Time, Telescope, TypedSet
from caom2.dali import Interval


def main():

    writer_25 = ow_25(validate=False, write_empty_collections=False)
    reader_24 = or_24(validate=False)
    obs_24 = reader_24.read(sys.argv[1])

    proposal = Proposal(
        id=obs_24.proposal.id,
        pi=obs_24.proposal.pi_name,
        title=obs_24.proposal.title,
    )
    for keyword in obs_24.proposal.keywords:
        proposal.keywords.add(keyword)
    target = Target(
        name=obs_24.target.name,
        type=obs_24.target.target_type.value,
    )
    telescope = Telescope(
        name=obs_24.telescope.name,
        geo_location_x=obs_24.telescope.geo_location_x,
        geo_location_y=obs_24.telescope.geo_location_y,
        geo_location_z=obs_24.telescope.geo_location_z,
    )
    instrument = Instrument(
        name=obs_24.instrument.name
    )
    obs = SimpleObservation(
        collection=obs_24.collection,
        uri=f'caom:{obs_24.collection}/{obs_24.observation_id}',
        algorithm=Algorithm(
            obs_24.algorithm.name
        ),
        proposal=proposal,
        target=target,
        telescope=telescope,
        instrument=instrument,
        intent=ObservationIntentType(obs_24.intent.value),
        type=obs_24.type,
        meta_release=obs_24.meta_release,
    )
    for plane_24 in obs_24.planes.values():
        provenance = Provenance(
            name=plane_24.provenance.name,
            version=plane_24.provenance.version,
            project=plane_24.provenance.project,
            producer=plane_24.provenance.producer,
            run_id=plane_24.provenance.run_id,
            reference=plane_24.provenance.reference,
        )
        for keyword in plane_24.provenance.keywords:
            provenance.keywords.add(keyword)
        for inpt in plane_24.provenance.inputs:
            provenance.inputs.add(inpt.uri)
        plane = Plane(
            uri=f'{obs.uri}/{plane_24.product_id}',
            meta_release=plane_24.meta_release,
            data_release=plane_24.data_release,
            data_product_type=DataProductType(plane_24.data_product_type.value),
            calibration_level=CalibrationLevel(plane_24.calibration_level.value),
            provenance=provenance,
        )
        # Erik Osinga
        # LoTSS data are taken (and stored) at a time resolution of 1 s and a frequency resolution of
        # 16 channels per 195.3 kHz sub-band (SB) by the observatory. Thus 12.1875 kHz frequency
        # resolution. The total bandwidth is 48 MHz, so that would be almost 4000 channels.
        # But for the LoTSS data (6'' resolution) we don't need such high resolution in time and
        # frequency. So  generally for imaging the LoTSS data people average to 8s time resolution
        # and 97.6 kHz frequency resolution. This is also the frequency resolution at which the
        # Stokes IQU cubes are produced. These cubes have 480 "channels" in that sense.
        num_channels = 480
        channel_width_in_mhz = 0.1  # MHz => ( 168 MHz - 120 MHz ) = 48 MHz => / 480 channels
        channel_width_in_m = (channel_width_in_mhz * units.Hz).to(units.m, equivalencies=units.spectral()).value
        min_exposure = None
        max_exposure = None
        for artifact_24 in plane_24.artifacts.values():
            for part_24 in artifact_24.parts.values():
                for chunk_24 in part_24.chunks:
                    if chunk_24.energy:
                        # units are 'm'
                        if min_exposure is None:
                            min_exposure = chunk_24.energy.axis.range.start.val
                        else:
                            min_exposure = min(min_exposure, chunk_24.energy.axis.range.start.val)
                        if max_exposure is None:
                            max_exposure = chunk_24.energy.axis.range.end.val
                        else:
                            max_exposure = max(max_exposure, chunk_24.energy.axis.range.end.val)
        resolution_bounds = None
        if min_exposure and max_exposure:
            if min_exposure > max_exposure:
                resolution_bounds = Interval(
                    # TB - energy axis resolution == channel width in CAOM units
                    # double-check with Pat that this is consistent with the intention of the model
                    # check with Pat that this value makes sense with LoTSS data and that resolution == channel width
                    # does the channel width vary as a function of frequency
                    lower=channel_width_in_m,
                    upper=channel_width_in_m,
                )
        energy_bounds = Interval(
            lower=plane_24.energy.bounds.lower,
            upper=plane_24.energy.bounds.upper,
        )
        energy_samples = []
        for energy_sample_24 in plane_24._energy.bounds.samples:
            energy_sample = Interval(
                lower=energy_sample_24.lower,
                upper= energy_sample_24.upper,
            )
            energy_samples.append(energy_sample)
        energy_bands = TypedSet(EnergyBand)
        for energy_band_24 in plane_24._energy.energy_bands:
            energy_bands.add(EnergyBand(energy_band_24.value))
        energy = Energy(
            bounds=energy_bounds,
            samples=energy_samples,
            dimension=plane_24._energy.dimension,
            sample_size=plane_24._energy.sample_size,
            energy_bands=energy_bands,
            rest=plane_24._energy.restwav,
            bandpass_name=plane_24._energy.bandpass_name,
            resolution_bounds=resolution_bounds,
        )
        plane._energy = energy
        if(
            plane_24._position
            and plane_24._position.bounds
            and isinstance(plane_24._position.bounds, Circle_24)
            and plane_24._position.bounds.center
        ):
            center_point = Point(
                cval1=plane_24._position.bounds.center.cval1,
                cval2=plane_24._position.bounds.center.cval2,
            )
            position_samples = Circle(
                center=center_point,
                radius=plane_24._position.bounds.radius,
            )
            position_bounds = position_samples

            position_dimension = Dimension2D(
                naxis1=plane_24._position.dimension.naxis1,
                naxis2=plane_24._position.dimension.naxis2,
            )
            # for maxRecoverableScale
            # https://www.aanda.org/articles/aa/pdf/2013/08/aa20873-12.pdf, Table 1
            # minimum baseline (m) 68
            # JJK
            # I think that for LoTSS they have done work to actually calibrate what the max
            # scale is.  From that work they find values of about 120'' (120/3600 degrees).
            # That analysis does not dwell on the frequency dependence as there are significant
            # noise and data processing issues that are more substantial than variation in
            # scale based on perfect systems.
            #
            # the model allows for a range in resolution intervals and a range in max recoverable scale. Do the
            # channel width and the mrs of the lotss data vary with frequency, in which case we can represent that
            # as an interval in the data model
            # NO => one value as lower, upper
            # YES => multiple values
            max_recoverable_scale = Interval(
                lower=120.0,
                upper=120.0,
            )
            position = Position(
                bounds=position_samples,
                min_bounds=position_samples,
                samples=MultiShape([position_samples]),
                dimension=position_dimension,
                sample_size=plane_24._position.sample_size,
                resolution=plane_24._position.resolution,
                max_recoverable_scale = max_recoverable_scale,
            )
            plane._position = position
        elif(
            plane_24._position
            and plane_24._position.bounds
            and isinstance(plane_24._position.bounds, Position_24)
        ):
            position_bounds = []
            for point_24 in plane_24._position.bounds.points:
                point = Point(
                    cval1=point_24.cval1,
                    cval2=point_24.cval2,
                )
                position_bounds.append(point)
            position_dimension = Dimension2D(
                naxis1=plane_24._position.dimension.naxis1,
                naxis2=plane_24._position.dimension.naxis2,
            )
            max_recoverable_scale = None
            position_polygon = Polygon(points=position_bounds)
            position_multi_shape = MultiShape(shapes=[position_polygon])
            position = Position(
                bounds=position_polygon,
                samples=position_multi_shape,
                min_bounds=position_polygon,
                dimension=position_dimension,
                max_recoverable_scale = max_recoverable_scale,
                resolution=plane_24._position.resolution,
                resolution_bounds=None,
                sample_size=plane_24._position.sample_size,
                calibration=None,
            )
            plane._position = position

        time_samples = None
        if plane_24._time and plane_24._time.bounds:
            time_samples = Interval(
                lower=plane_24._time.bounds.lower,
                upper=plane_24._time.bounds.upper,
            )

            min_duration = None
            max_duration = None
            for artifact_24 in plane_24.artifacts.values():
                for part_24 in artifact_24.parts.values():
                    for chunk_24 in part_24.chunks:
                        if chunk_24.time and chunk_24.time.exposure:
                            if min_duration is None:
                                min_duration = chunk_24.time.exposure
                            else:
                                min_duration = min(min_duration, chunk_24.time.exposure)
                            if max_duration is None:
                                max_duration = chunk_24.time.exposure
                            else:
                                max_duration = max(max_duration, chunk_24.time.exposure)
            exposure_bounds = None
            if min_duration and max_duration:
                exposure_bounds = Interval(
                    lower=min_duration,
                    upper=max_duration,
                )
            time = Time(
                bounds=time_samples,
                samples=[time_samples],
                dimension=plane_24._time.dimension,
                sample_size=plane_24._time.sample_size,
                resolution=plane_24._time.resolution,
                exposure=plane_24._time.exposure,
                exposure_bounds=exposure_bounds,
            )
            plane._time = time

        if plane_24._polarization and plane_24._polarization.dimension and plane_24._polarization.polarization_states:
            pol_states = []
            for state in plane_24._polarization.polarization_states:
                pol_states.append(PolarizationState(state.value))
            polarization = Polarization(
                dimension=plane_24._polarization.dimension,
                states=plane_24._polarization.polarization_states,
            )
            plane._polarization = polarization
        obs.planes.add(plane)

        for artifact_24 in plane_24.artifacts.values():
            artifact = Artifact(
                uri=artifact_24.uri,
                product_type=DataLinkSemantics(artifact_24.product_type.value),
                release_type=ReleaseType(artifact_24.release_type.value),
                content_type=artifact_24.content_type,
                content_length=artifact_24.content_length,
                content_checksum=artifact_24.content_checksum.uri if artifact_24.content_checksum else None,
            )
            plane.artifacts.add(artifact)

    output_f_name = sys.argv[1].replace('.24.', '.25.').replace('inputs', 'outputs')
    writer_25.write(obs, output_f_name)




if __name__ == '__main__':
    try:
        print('start')
        main()
        print('done')
    except Exception as e:
        import traceback
        print(traceback.format_exc())
        print(e)

