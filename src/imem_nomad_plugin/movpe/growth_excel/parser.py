#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

from time import sleep, perf_counter
import pandas as pd
import yaml
import json
from typing import Dict, List

from nomad.units import ureg
from nomad.utils import hash
from nomad.metainfo import MSection, Quantity, Section
from nomad.parsing import MatchingParser
from nomad.datamodel.metainfo.annotations import ELNAnnotation
from nomad.datamodel.data import EntryData
from nomad.datamodel.datamodel import EntryArchive, EntryMetadata
from nomad.datamodel.metainfo.basesections import (
    SystemComponent,
    PubChemPureSubstanceSection,
    CompositeSystemReference,
)

from nomad_material_processing import (
    SubstrateReference,
    ThinFilmReference,
    ThinFilmStackReference,
    Parallelepiped,
    SubstrateCrystalProperties,
    Miscut,
    Dopant,
)
from nomad_material_processing.vapor_deposition import Pressure, VolumetricFlowRate
from nomad_material_processing.vapor_deposition.cvd import Rotation

from imem_nomad_plugin.general.schema import (
    SampleCutIMEM,
)

from imem_nomad_plugin.characterization.schema import (
    AFMmeasurement,
    AFMresults,
)

from imem_nomad_plugin.movpe.schema import (
    ExperimentMovpeIMEM,
    GrowthStepMovpeIMEM,
    GrowthMovpeIMEM,
    GrowthMovpeIMEMReference,
    SampleCutIMEMReference,
    ThinFilmMovpe,
    ThinFilmStackMovpeIMEM,
    ThinFilmStackMovpeReference,
    SampleParametersMovpe,
    ChamberEnvironmentMovpe,
    GasFlowMovpe,
    ShaftTemperature,
    FilamentTemperature,
    SubstrateMovpe,
    SubstrateCrystalPropertiesMovpe,
    MiscutMovpe,
    Shape,
    SampleParametersMovpe,
    FilamentTemperature,
    XRDmeasurementReference,
    AFMmeasurementReference,
    CharacterizationMovpeIMEM,
)

from imem_nomad_plugin.utils import (
    create_archive,
    fetch_substrate,
    populate_sources,
    populate_gas_source,
    fill_quantity,
    populate_element,
    populate_dopant,
)


class RawFileGrowthRun(EntryData):
    m_def = Section(a_eln=None, label='Raw File Growth Run')
    name = Quantity(
        type=str,
        a_eln=ELNAnnotation(
            component='StringEditQuantity',
        ),
    )
    growth_run = Quantity(
        type=ExperimentMovpeIMEM,
        # a_eln=ELNAnnotation(
        #     component="ReferenceEditQuantity",
        # ),
        shape=['*'],
    )


class ParserMovpeIMEM(MatchingParser):
    def parse(self, mainfile: str, archive: EntryArchive, logger) -> None:
        xlsx = pd.ExcelFile(mainfile)
        data_file = mainfile.split('/')[-1]
        data_file_with_path = mainfile.split('raw/')[-1]
        filetype = 'yaml'

        # creates a new DataFrame with the stripped column names.
        # sheet = pd.read_excel(xlsx, 'Overview', comment='#', converters={'Sample': str})
        # overview_sheet = sheet.rename(columns=lambda x: x.strip())

        if 'Overview' in xlsx.sheet_names:
            overview_sheet = pd.read_excel(
                xlsx, 'Overview', comment='#', converters={'Sample': str}
            )
            overview_sheet.columns = overview_sheet.columns.str.strip()
        if 'Substrate' in xlsx.sheet_names:
            substrates_sheet = pd.read_excel(
                xlsx,
                'Substrate',
                comment='#',
                converters={'Orientation': str, 'Off-cut Orientation': str},
            )
            substrates_sheet.columns = substrates_sheet.columns.str.strip()
        if 'GrowthRun' in xlsx.sheet_names:
            growthrun_sheet = pd.read_excel(xlsx, 'GrowthRun', comment='#')
            growthrun_sheet.columns = growthrun_sheet.columns.str.strip()
        if 'Precursors' in xlsx.sheet_names:
            precursors_sheet = pd.read_excel(xlsx, 'Precursors', comment='#')
            precursors_sheet.columns = precursors_sheet.columns.str.strip()
        if 'Mist' in xlsx.sheet_names:
            mist_sheet = pd.read_excel(xlsx, 'Mist', comment='#')
            mist_sheet.columns = mist_sheet.columns.str.strip()
        if 'Pregrowth' in xlsx.sheet_names:
            pregrowth_sheet = pd.read_excel(xlsx, 'Pregrowth', comment='#')
            pregrowth_sheet.columns = pregrowth_sheet.columns.str.strip()
        if 'SampleCut' in xlsx.sheet_names:
            samplecut_sheet = pd.read_excel(xlsx, 'SampleCut', comment='#')
            samplecut_sheet.columns = samplecut_sheet.columns.str.strip()
        if 'HRXRD' in xlsx.sheet_names:
            hrxrd_sheet = pd.read_excel(
                xlsx, 'HRXRD', converters={'Sample': str}, comment='#'
            )
            hrxrd_sheet.columns = hrxrd_sheet.columns.str.strip()
        if 'AFMReflectanceSEM' in xlsx.sheet_names:
            characterization_sheet = pd.read_excel(
                xlsx, 'AFMReflectanceSEM', converters={'Sample': str}, comment='#'
            )
            characterization_sheet.columns = characterization_sheet.columns.str.strip()
        if 'ElectroOptical' in xlsx.sheet_names:
            electro_optical_sheet = pd.read_excel(xlsx, 'ElectroOptical', comment='#')
            electro_optical_sheet.columns = electro_optical_sheet.columns.str.strip()
        if 'Contacts' in xlsx.sheet_names:
            contacts_sheet = pd.read_excel(xlsx, 'Contacts', comment='#')
            contacts_sheet.columns = contacts_sheet.columns.str.strip()

        try:
            sample_id = fill_quantity(overview_sheet.iloc[0], 'Sample')
        except IndexError:
            sample_id = None

        # creating Substrate archives
        substrate_index = None
        for (
            substrate_index,
            substrate_row,
        ) in substrates_sheet.iterrows():
            substrate_id = fill_quantity(substrate_row, 'Substrates')
            # creating Substrate archives
            substrate_filename = (
                f'{substrate_id}_{substrate_index}.SubstrateIKZ.archive.{filetype}'
            )
            substrate_data = SubstrateMovpe()
            substrate_data.geometry = Shape()
            substrate_data.crystal_properties = SubstrateCrystalPropertiesMovpe()
            substrate_data.crystal_properties.miscut = MiscutMovpe()

            substrate_data.lab_id = substrate_id
            substrate_data.name = fill_quantity(substrate_row, 'Material')

            de = fill_quantity(substrate_row, 'Description')
            notes = fill_quantity(substrate_row, 'Notes')
            substrate_data.description = f'Description: {de}, Notes: {notes}'

            substrate_supplier_id = fill_quantity(substrate_row, 'Substrate ID')
            if substrate_supplier_id:
                substrate_data.supplier_id = substrate_supplier_id

            substrate_supplier = fill_quantity(substrate_row, 'Supplier')
            if substrate_supplier:
                substrate_data.supplier = substrate_supplier

            substrate_data_geometry_width = fill_quantity(
                substrate_row,
                'Size X',
                read_unit='millimeter',
            )
            if substrate_data_geometry_width:
                substrate_data.geometry.width = substrate_data_geometry_width

            substrate_data_geometry_length = fill_quantity(
                substrate_row,
                'Size Y',
                read_unit='millimeter',
            )
            if substrate_data_geometry_length:
                substrate_data.geometry.length = substrate_data_geometry_length

            substrate_data_geometry_diameter = fill_quantity(
                substrate_row,
                'Size Diameter',
                read_unit='millimeter',
            )
            if substrate_data_geometry_diameter:
                substrate_data.geometry.diameter = substrate_data_geometry_diameter

            annealing = fill_quantity(substrate_row, 'Annealing')
            if annealing:
                substrate_data.annealing = annealing

            cleaning = fill_quantity(substrate_row, 'Cleaning')
            if cleaning:
                substrate_data.cleaning = cleaning

            regrowth = fill_quantity(substrate_row, 'Regrowth')
            if regrowth:
                substrate_data.regrowth = regrowth

            orientation = fill_quantity(substrate_row, 'Orientation')
            if orientation:
                substrate_data.crystal_properties.orientation = orientation

            miscut_angle = fill_quantity(substrate_row, 'Off-cut')
            if miscut_angle:
                substrate_data.crystal_properties.miscut.angle = miscut_angle

            mo = fill_quantity(substrate_row, 'Off-cut Orientation')
            if mo:
                substrate_data.crystal_properties.miscut.orientation = mo

            # TODO check functions populate_element and populate_dopant
            # and make them compliant to the other parsing
            substrate_data.elemental_composition = populate_element(
                substrate_index, substrates_sheet
            )

            substrate_data.dopants = populate_dopant(substrate_index, substrates_sheet)

            substrate_archive = EntryArchive(
                data=substrate_data if substrate_data else SubstrateMovpe(),
                m_context=archive.m_context,
                metadata=EntryMetadata(upload_id=archive.m_context.upload_id),
            )
            create_archive(
                substrate_archive.m_to_dict(),
                archive.m_context,
                substrate_filename,
                filetype,
                logger,
            )

        # generate substrate references
        try:
            substrate_id = fill_quantity(substrates_sheet.iloc[0], 'Substrates')
        except IndexError:
            substrate_id = None
        # TODO try to get rid of the fetch_substarte as it slows down the processing
        sleep(1.5)
        substrate_ref = None
        if substrate_id and sample_id:
            substrate_ref = fetch_substrate(archive, sample_id, substrate_id, logger)
        sub_ref = None
        if substrate_ref is not None:
            sub_ref = SubstrateReference(reference=substrate_ref)
        else:
            sub_ref = SubstrateReference(name=substrate_id, lab_id=substrate_id)

        # creating ThinFiln and ThinFilmStack archives
        layer_filename = f'{sample_id}.ThinFilm.archive.{filetype}'
        layer_archive = EntryArchive(
            data=ThinFilmMovpe(
                name=str(sample_id) + ' layer',
                lab_id=str(sample_id) + 'layer',
            ),
            m_context=archive.m_context,
            metadata=EntryMetadata(upload_id=archive.m_context.upload_id),
        )
        create_archive(
            layer_archive.m_to_dict(),
            archive.m_context,
            layer_filename,
            filetype,
            logger,
        )
        grown_sample_data = ThinFilmStackMovpeIMEM(
            name=str(sample_id) + ' stack',
            lab_id=str(sample_id),
            layers=[
                ThinFilmReference(
                    reference=f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.m_context.upload_id, layer_filename)}#data'
                )
            ],
            substrate=sub_ref,
        )

        grown_sample_filename = f'{sample_id}.ThinFilmStack.archive.{filetype}'
        grown_sample_archive = EntryArchive(
            data=grown_sample_data if grown_sample_data else ThinFilmStackMovpeIMEM(),
            m_context=archive.m_context,
            metadata=EntryMetadata(upload_id=archive.m_context.upload_id),
        )
        create_archive(
            grown_sample_archive.m_to_dict(),
            archive.m_context,
            grown_sample_filename,
            filetype,
            logger,
        )

        # creating samplecut archive
        samplecut_data = SampleCutIMEM(
            name=str(sample_id) + ' sample cut',
            parent_sample=ThinFilmStackMovpeReference(
                reference=f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.m_context.upload_id, grown_sample_filename)}#data'
            ),
        )
        for _, row in samplecut_sheet.iterrows():
            children_sample_id = (
                f"{sample_id} {fill_quantity(row, 'Children Sample ID')}"
            )
            children_sample_data = ThinFilmStackMovpeIMEM(
                name=str(children_sample_id) + ' stack',
                lab_id=str(children_sample_id),
                layers=[
                    ThinFilmReference(
                        reference=f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.m_context.upload_id, layer_filename)}#data'
                    )
                ],
                substrate=sub_ref,
                parent_sample=ThinFilmStackMovpeReference(
                    reference=f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.m_context.upload_id, grown_sample_filename)}#data'
                ),
            )
            children_sample_filename = (
                f'{children_sample_id}.ThinFilmStack.archive.{filetype}'
            )
            children_sample_archive = EntryArchive(
                data=children_sample_data
                if children_sample_data
                else ThinFilmStackMovpeIMEM(),
                m_context=archive.m_context,
                metadata=EntryMetadata(upload_id=archive.m_context.upload_id),
            )
            create_archive(
                children_sample_archive.m_to_dict(),
                archive.m_context,
                children_sample_filename,
                filetype,
                logger,
            )
            samplecut_data.children_samples.append(
                ThinFilmStackMovpeReference(
                    reference=f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.m_context.upload_id, children_sample_filename)}#data'
                )
            )
        samplecut_filename = f'{sample_id}.SampleCut.archive.{filetype}'
        samplecut_archive = EntryArchive(
            data=samplecut_data if samplecut_data else SampleCutIMEM(),
            m_context=archive.m_context,
            metadata=EntryMetadata(upload_id=archive.m_context.upload_id),
        )
        create_archive(
            samplecut_archive.m_to_dict(),
            archive.m_context,
            samplecut_filename,
            filetype,
            logger,
        )

        # creating growth process step objects
        process_steps_lists = []
        for step_index, step in growthrun_sheet.iterrows():
            growth_step = GrowthStepMovpeIMEM()
            growth_step.environment = ChamberEnvironmentMovpe()
            growth_step.environment.pressure = Pressure()
            growth_step.environment.rotation = Rotation()
            growth_step.environment.gas_flow = [GasFlowMovpe()]
            growth_step.environment.gas_flow[0].gas = PubChemPureSubstanceSection()
            growth_step.environment.uniform_gas_flow_rate = VolumetricFlowRate()
            growth_step.sample_parameters = [SampleParametersMovpe()]
            growth_step.sample_parameters[
                0
            ].filament_temperature = FilamentTemperature()

            growth_step.step_id = step_index + 1
            growth_step.comment = fill_quantity(step, 'Notes')

            growth_step_name = fill_quantity(step, 'Name')
            if growth_step_name:
                growth_step.name = growth_step_name

            growth_step_duration = fill_quantity(step, 'Duration', read_unit='minute')
            if growth_step_duration:
                growth_step.duration = growth_step_duration

            growth_step_environment_pressure_set_value = fill_quantity(
                step, 'Pressure', read_unit='mbar'
            )
            if growth_step_environment_pressure_set_value:
                growth_step.environment.pressure.set_time = pd.Series([0])
                growth_step.environment.pressure.set_value = pd.Series(
                    [growth_step_environment_pressure_set_value]
                )

            growth_step_environment_rotation_set_value = fill_quantity(
                step, 'Rotation', read_unit='rpm'
            )
            if growth_step_environment_rotation_set_value:
                growth_step.environment.rotation.set_time = pd.Series([0])
                growth_step.environment.rotation.set_value = pd.Series(
                    [growth_step_environment_rotation_set_value]
                )

            growth_step_environment_gas_flow_gas_name = fill_quantity(
                step, 'Carrier Gas'
            )
            if growth_step_environment_gas_flow_gas_name:
                growth_step.environment.gas_flow[
                    0
                ].gas.name = growth_step_environment_gas_flow_gas_name

            growth_step_environment_uniform_gas_flow_rate_set_value = fill_quantity(
                step, 'Uniform Valve', read_unit='cm ** 3 / minute'
            )
            if growth_step_environment_uniform_gas_flow_rate_set_value:
                growth_step.environment.uniform_gas_flow_rate.set_time = pd.Series([0])
                growth_step.environment.uniform_gas_flow_rate.set_value = pd.Series(
                    [growth_step_environment_uniform_gas_flow_rate_set_value]
                )

            growth_step_sample_parameters_filament_temperature_set_value = (
                fill_quantity(step, 'Temperature', read_unit='celsius')
            )
            if growth_step_sample_parameters_filament_temperature_set_value:
                growth_step.sample_parameters[
                    0
                ].filament_temperature.set_time = pd.Series([0])
                growth_step.sample_parameters[
                    0
                ].filament_temperature.set_value = pd.Series(
                    [growth_step_sample_parameters_filament_temperature_set_value]
                )

            growth_step.sources = populate_sources(
                step_index, growthrun_sheet
            ) + populate_gas_source(step_index, growthrun_sheet)

            process_steps_lists.append(growth_step)

        # creating growth process objects
        growth_process_object = GrowthMovpeIMEM()

        growth_process_object.lab_id = sample_id

        if not substrates_sheet.empty:
            growth_susceptor = fill_quantity(substrates_sheet.iloc[0], 'Susceptor')
            if growth_susceptor:
                growth_process_object.susceptor = growth_susceptor

            growth_mask = fill_quantity(substrates_sheet.iloc[0], 'Mask')
            if growth_mask:
                growth_process_object.mask = growth_mask

            growth_pocket = fill_quantity(substrates_sheet.iloc[0], 'Pocket')
            if growth_pocket:
                growth_process_object.pocket = growth_pocket

        growth_process_object.samples = [
            ThinFilmStackMovpeReference(
                reference=f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.m_context.upload_id, grown_sample_filename)}#data'
            ),
        ]

        growth_process_object.steps = process_steps_lists

        # creating growth process archives
        growth_process_filename = f'{sample_id}.GrowthMovpeIMEM.archive.{filetype}'
        # Activity.normalize(growth_process_object, archive, logger)
        growth_process_archive = EntryArchive(
            data=growth_process_object if growth_process_object else GrowthMovpeIMEM(),
            m_context=archive.m_context,
            metadata=EntryMetadata(upload_id=archive.m_context.upload_id),
        )

        create_archive(
            growth_process_archive.m_to_dict(),
            archive.m_context,
            growth_process_filename,
            filetype,
            logger,
        )

        # creating growth PRE-process step objects
        pre_process_steps_lists = []
        for step_id, step in pregrowth_sheet.iterrows():
            pregrowth_step = GrowthStepMovpeIMEM()
            pregrowth_step.environment = ChamberEnvironmentMovpe()
            pregrowth_step.environment.pressure = Pressure()
            pregrowth_step.environment.rotation = Rotation()
            pregrowth_step.environment.gas_flow = [GasFlowMovpe()]
            pregrowth_step.environment.gas_flow[0].gas = PubChemPureSubstanceSection()
            pregrowth_step.environment.uniform_gas_flow_rate = VolumetricFlowRate()
            pregrowth_step.sample_parameters = [SampleParametersMovpe()]
            pregrowth_step.sample_parameters[
                0
            ].filament_temperature = FilamentTemperature()

            pregrowth_step.step_id = step_id + 1
            pregrowth_step.comment = fill_quantity(step, 'Notes')

            pregrowth_step_name = fill_quantity(step, 'Step Name')
            if pregrowth_step_name:
                pregrowth_step.name = pregrowth_step_name

            pregrowth_step_duration = fill_quantity(
                step, 'Duration', read_unit='minute'
            )
            if pregrowth_step_duration:
                pregrowth_step.duration = pregrowth_step_duration

            pregrowth_step_environment_pressure_set_value = fill_quantity(
                step, 'Chamber Pressure', read_unit='mbar'
            )
            if pregrowth_step_environment_pressure_set_value:
                pregrowth_step.environment.pressure.set_time = pd.Series([0])
                pregrowth_step.environment.pressure.set_value = pd.Series(
                    [pregrowth_step_environment_pressure_set_value]
                )

            pregrowth_step_environment_rotation_set_value = fill_quantity(
                step, 'Carrier Rotation', read_unit='rpm'
            )
            if pregrowth_step_environment_rotation_set_value:
                pregrowth_step.environment.rotation.set_time = pd.Series([0])
                pregrowth_step.environment.rotation.set_value = pd.Series(
                    [pregrowth_step_environment_rotation_set_value]
                )

            pregrowth_step_environment_gas_flow_gas_name = fill_quantity(
                step, 'Carrier Gas'
            )
            if pregrowth_step_environment_gas_flow_gas_name:
                pregrowth_step.environment.gas_flow[
                    0
                ].gas.name = pregrowth_step_environment_gas_flow_gas_name

            pregrowth_step_environment_uniform_gas_flow_rate_set_value = fill_quantity(
                step, 'Carrier Gas Flow', read_unit='cm ** 3 / minute'
            )
            if pregrowth_step_environment_uniform_gas_flow_rate_set_value:
                pregrowth_step.environment.uniform_gas_flow_rate.set_time = pd.Series(
                    [0]
                )
                pregrowth_step.environment.uniform_gas_flow_rate.set_value = pd.Series(
                    [pregrowth_step_environment_uniform_gas_flow_rate_set_value]
                )

            pregrowth_step_sample_parameters_filament_temperature_set_value = (
                fill_quantity(step, 'Substrate Temperature', read_unit='celsius')
            )
            if pregrowth_step_sample_parameters_filament_temperature_set_value:
                pregrowth_step.sample_parameters[
                    0
                ].filament_temperature.set_time = pd.Series([0])
                pregrowth_step.sample_parameters[
                    0
                ].filament_temperature.set_value = pd.Series(
                    [pregrowth_step_sample_parameters_filament_temperature_set_value]
                )
            pre_process_steps_lists.append(pregrowth_step)

        # creating PRE-growth process objects
        pregrowth_process_object = GrowthMovpeIMEM()

        if not substrates_sheet.empty:
            growth_susceptor = fill_quantity(substrates_sheet.iloc[0], 'Susceptor')
            if not growth_susceptor:
                pregrowth_process_object.susceptor = growth_susceptor

            growth_mask = fill_quantity(substrates_sheet.iloc[0], 'Mask')
            if not growth_mask:
                pregrowth_process_object.mask = growth_mask

            growth_pocket = fill_quantity(substrates_sheet.iloc[0], 'Pocket')
            if not growth_pocket:
                pregrowth_process_object.pocket = growth_pocket

        pregrowth_process_object.lab_id = sample_id

        # creating pregrowth process archives
        pregrowth_process_filename = (
            f'{sample_id}-pregrowth.GrowthMovpeIMEM.archive.{filetype}'
        )
        # Activity.normalize(pregrowth_process_object, archive, logger)
        pregrowth_process_archive = EntryArchive(
            data=pregrowth_process_object
            if pregrowth_process_object
            else GrowthMovpeIMEM(),
            m_context=archive.m_context,
            metadata=EntryMetadata(upload_id=archive.m_context.upload_id),
        )

        create_archive(
            pregrowth_process_archive.m_to_dict(),
            archive.m_context,
            pregrowth_process_filename,
            filetype,
            logger,
        )

        # CONTINUA DA QUA

        # sheet = pd.read_excel(
        #     xlsx, 'ElectroOptical', comment='#', converters={'Sample': str}
        # )
        # hall_measurement_sheet = sheet.rename(columns=lambda x: x.strip())
        # hall_meas_refs = []
        # for meas_index, hall_measurement in hall_measurement_sheet.iterrows():
        #     filetype = 'yaml'
        #     hall_meas_filename = (
        #         f"{hall_measurement['Sample']}_{meas_index}_hall.archive.{filetype}"
        #     )
        #     hall_meas_archive = EntryArchive(
        #         data=HallMeasurement(
        #             lab_id=hall_measurement['Sample'],
        #             datetime=datetime.datetime.strptime(
        #                 hall_measurement['Date'], r'%Y-%m-%d'
        #             ).astimezone(),
        #             results=[
        #                 HallMeasurementResult(
        #                     resistivity=hall_measurement['Resistivity'],
        #                     mobility=hall_measurement['Mobility'],
        #                     carrier_concentration=hall_measurement['Carrier Concentration'],
        #                 )
        #             ],
        #         ),
        #         m_context=archive.m_context,
        #         metadata=EntryMetadata(upload_id=archive.m_context.upload_id),
        #     )
        #     create_archive(
        #         hall_meas_archive.m_to_dict(),
        #         archive.m_context,
        #         hall_meas_filename,
        #         filetype,
        #         logger,
        #     )
        #     hall_meas_refs.append(
        #         HallMeasurements(
        #             reference=f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.metadata.upload_id, hall_meas_filename)}#data'
        #         )
        #     )

        # creating AFM archive
        afm_measurements = []
        for index, row in characterization_sheet.iterrows():
            afm_data = AFMmeasurement()
            afm_data.samples = []
            afm_data.results = [AFMresults()]

            afm_name = fill_quantity(row, 'Sample')
            if afm_name:
                afm_data.name = f'{afm_name} afm {index}'
                afm_data.samples.append(
                    CompositeSystemReference(
                        lab_id=str(afm_name),
                    )
                )
                afm_filename = f'{afm_name}_{index}.AFMmeasurement.archive.{filetype}'
            else:
                afm_filename = f'{sample_id}_{index}.AFMmeasurement.archive.{filetype}'

            afm_datetime = fill_quantity(row, 'Date')
            if afm_datetime:
                afm_data.datetime = afm_datetime

            afm_results_notes = fill_quantity(row, 'Notes')
            if afm_results_notes:
                afm_data.results[0].name = afm_results_notes

            afm_results_roughness = fill_quantity(row, 'Roughness')
            if afm_results_roughness:
                afm_data.results[0].roughness = afm_results_roughness

            afm_results_surface_features = fill_quantity(row, 'Surface Features')
            if afm_results_surface_features:
                afm_data.results[0].surface_features = afm_results_surface_features

            afm_archive = EntryArchive(
                data=afm_data if afm_data else AFMmeasurement(),
                m_context=archive.m_context,
                metadata=EntryMetadata(upload_id=archive.m_context.upload_id),
            )
            create_archive(
                afm_archive.m_to_dict(),
                archive.m_context,
                afm_filename,
                filetype,
                logger,
            )
            afm_measurements.append(
                AFMmeasurementReference(
                    name=f'{afm_name} afm {index}',
                    reference=f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.m_context.upload_id, afm_filename)}#data',
                )
            )

        # creating XRD object(s)
        xrd_measurements = []
        for index, row in hrxrd_sheet.iterrows():
            xrd_measurement = XRDmeasurementReference()

            xrd_name = fill_quantity(row, 'Sample')
            if xrd_name:
                xrd_measurement.name = f'{xrd_name} xrd {index}'

            xrd_phase = fill_quantity(row, 'Phase')
            if xrd_phase:
                xrd_measurement.phase = xrd_phase

            xrd_peak_position_2theta = fill_quantity(row, 'Peak Position - 2theta')
            if xrd_peak_position_2theta:
                xrd_measurement.peak_position_2theta = xrd_peak_position_2theta

            xrd_peak_fwhm_2theta = fill_quantity(row, 'Peak FWHM - 2theta')
            if xrd_peak_fwhm_2theta:
                xrd_measurement.peak_fwhm_2theta = xrd_peak_fwhm_2theta

            xrd_peak_position_omega = fill_quantity(row, 'Peak Position - Omega')
            if xrd_peak_position_omega:
                xrd_measurement.peak_position_omega = xrd_peak_position_omega

            xrd_peak_fwhm_rocking_curve = fill_quantity(row, 'Peak FWHM Rocking Curve')
            if xrd_peak_fwhm_rocking_curve:
                xrd_measurement.peak_fwhm_rocking_curve = xrd_peak_fwhm_rocking_curve

            xrd_reflection = fill_quantity(row, 'Reflection')
            if xrd_reflection:
                xrd_measurement.reflection = xrd_reflection

            xrd_description = fill_quantity(row, 'Notes')
            if xrd_description:
                xrd_measurement.description = xrd_description

            xrd_measurements.append(xrd_measurement)

        # creating experiment archive
        experiment_filename = f'{sample_id}.ExperimentMovpeIMEM.archive.{filetype}'
        experiment_data = ExperimentMovpeIMEM(
            name='experiment',
            method='MOVPE 2 experiment',
            pregrowth=GrowthMovpeIMEMReference(
                reference=f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.m_context.upload_id, pregrowth_process_filename)}#data',
            ),
            growth_run=GrowthMovpeIMEMReference(
                reference=f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.m_context.upload_id, growth_process_filename)}#data',
            ),
            sample_cut=SampleCutIMEMReference(
                reference=f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.m_context.upload_id, samplecut_filename)}#data',
            ),
            characterization=CharacterizationMovpeIMEM(
                xrd=xrd_measurements,
                afm=afm_measurements,
            ),
        )
        experiment_archive = EntryArchive(
            data=experiment_data if experiment_data else ExperimentMovpeIMEM(),
            # m_context=archive.m_context,
            metadata=EntryMetadata(upload_id=archive.m_context.upload_id),
        )
        create_archive(
            experiment_archive.m_to_dict(),
            archive.m_context,
            experiment_filename,
            filetype,
            logger,
        )

        archive.data = RawFileGrowthRun(
            name=data_file,
            growth_run=[
                f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.m_context.upload_id, experiment_filename)}#data'
            ],
        )
        archive.metadata.entry_name = data_file + 'raw file'
