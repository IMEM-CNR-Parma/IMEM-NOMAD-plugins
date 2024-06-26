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
    typed_df_value,
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

        sample_id = overview_sheet['Sample'][0]

        # creating Substrate archives
        for substrate_index, substrate_id in enumerate(substrates_sheet['Substrates']):
            # creating Substrate archives
            substrate_filename = (
                f'{substrate_id}_{substrate_index}.SubstrateIKZ.archive.{filetype}'
            )
            substrate_data = SubstrateMovpe(
                lab_id=substrate_id,
                supplier_id=typed_df_value(
                    substrates_sheet, 'Substrate ID', str, substrate_index
                ),
                supplier=(
                    typed_df_value(substrates_sheet, 'Supplier', str, substrate_index)
                ),
                name=typed_df_value(substrates_sheet, 'Material', str, substrate_index),
                description=f"Description: {typed_df_value(substrates_sheet, 'Description', str, substrate_index)}, Notes: {typed_df_value(substrates_sheet, 'Notes', str, substrate_index)}",
                geometry=Shape(
                    width=typed_df_value(
                        substrates_sheet, 'Size X', float, substrate_index
                    )
                    * ureg('millimeter').to('meter').magnitude,
                    length=typed_df_value(
                        substrates_sheet, 'Size Y', float, substrate_index
                    )
                    * ureg('millimeter').to('meter').magnitude,
                    diameter=typed_df_value(
                        substrates_sheet, 'Size Diameter', float, substrate_index
                    )
                    * ureg('millimeter').to('meter').magnitude,
                ),
                crystal_properties=SubstrateCrystalPropertiesMovpe(
                    orientation=(
                        typed_df_value(
                            substrates_sheet, 'Orientation', str, substrate_index
                        )
                    ),
                    miscut=MiscutMovpe(
                        angle=(
                            typed_df_value(
                                substrates_sheet, 'Off-cut', float, substrate_index
                            )
                        ),
                        orientation=(
                            typed_df_value(
                                substrates_sheet,
                                'Off-cut Orientation',
                                str,
                                substrate_index,
                            )
                        ),
                    ),
                ),
                elemental_composition=populate_element(
                    substrate_index, substrates_sheet
                ),
                dopants=populate_dopant(substrate_index, substrates_sheet),
                annealing=typed_df_value(
                    substrates_sheet, 'Annealing', str, substrate_index
                ),
                cleaning=typed_df_value(
                    substrates_sheet, 'Cleaning', str, substrate_index
                ),
                regrowth=typed_df_value(
                    substrates_sheet, 'Regrowth', str, substrate_index
                ),
            )

            substrate_archive = EntryArchive(
                data=substrate_data,
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
        substrate_id = substrates_sheet['Substrates'][0]
        sleep(1.5)
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
            data=grown_sample_data,
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
            children_sample_id = row['Children Sample ID']
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
                data=children_sample_data,
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
            data=samplecut_data,
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
            process_steps_lists.append(
                GrowthStepMovpeIMEM(
                    name=(step['Name'] if 'Name' in growthrun_sheet.columns else None),
                    step_index=step_index + 1,
                    duration=(
                        step['Duration']
                        if 'Duration' in growthrun_sheet.columns
                        and not pd.isna(step['Duration'])
                        else 0 * ureg('minute').to('second').magnitude
                    ),
                    comment=(
                        step['Notes'] if 'Notes' in growthrun_sheet.columns else None
                    ),
                    sources=populate_sources(step_index, growthrun_sheet)
                    + populate_gas_source(step_index, growthrun_sheet),
                    environment=ChamberEnvironmentMovpe(
                        pressure=Pressure(
                            set_time=(pd.Series([0])),
                            set_value=pd.Series(
                                [
                                    (
                                        step['Pressure']
                                        if 'Pressure' in growthrun_sheet.columns
                                        else None
                                    )
                                ]
                            )
                            * ureg('mbar').to('pascal').magnitude,
                        ),
                        rotation=Rotation(
                            set_time=(pd.Series([0])),
                            set_value=pd.Series(
                                [
                                    (
                                        step['Rotation']
                                        if 'Rotation' in growthrun_sheet.columns
                                        else None
                                    )
                                ]
                            )
                            * ureg('rpm').to('rpm').magnitude,
                        ),
                        gas_flow=[
                            GasFlowMovpe(
                                gas=PubChemPureSubstanceSection(
                                    name=(
                                        step['Carrier Gas']
                                        if 'Carrier Gas' in growthrun_sheet.columns
                                        else None
                                    ),
                                ),
                            ),
                        ],
                        uniform_gas_flow_rate=VolumetricFlowRate(
                            set_time=(pd.Series([0])),
                            set_value=pd.Series(
                                [
                                    (
                                        step['Uniform Valve']
                                        if 'Uniform Valve' in growthrun_sheet.columns
                                        else None
                                    )
                                ]
                            )
                            * ureg('cm ** 3 / minute')
                            .to('meter ** 3 / second')
                            .magnitude,
                        ),
                    ),
                    sample_parameters=[
                        SampleParametersMovpe(
                            filament_temperature=FilamentTemperature(
                                time=(pd.Series([0])),
                                value=(
                                    pd.Series(
                                        [
                                            step['Temperature']
                                            if 'Temperature' in growthrun_sheet.columns
                                            else None
                                        ]
                                    )
                                    * ureg('celsius').to('kelvin').magnitude
                                ),
                                set_time=(pd.Series([0])),
                                set_value=(
                                    pd.Series(
                                        [
                                            step['Temperature']
                                            if 'Temperature' in growthrun_sheet.columns
                                            else None
                                        ]
                                    )
                                    * ureg('celsius').to('kelvin').magnitude
                                ),
                            )
                        )
                    ],
                )
            )

            # creating growth process objects
            growth_process_object = GrowthMovpeIMEM(
                name='Growth MOVPE',
                lab_id=sample_id,
                susceptor=typed_df_value(
                    substrates_sheet, 'Susceptor', str, substrate_index
                ),
                mask=typed_df_value(substrates_sheet, 'Mask', str, substrate_index),
                pocket=typed_df_value(substrates_sheet, 'Pocket', str, substrate_index),
                samples=[
                    ThinFilmStackMovpeReference(
                        reference=f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.m_context.upload_id, grown_sample_filename)}#data'
                    ),
                ],
            )
            growth_process_object.steps = process_steps_lists

        # creating growth process archives
        growth_process_filename = f'{sample_id}.GrowthMovpeIMEM.archive.{filetype}'
        # Activity.normalize(growth_process_object, archive, logger)
        growth_process_archive = EntryArchive(
            data=growth_process_object,
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
            pre_process_steps_lists.append(
                GrowthStepMovpeIMEM(
                    name=(
                        step['Step Name']
                        if 'Step Name' in pregrowth_sheet.columns
                        else None
                    ),
                    step_index=step_id + 1,
                    duration=(
                        step['Duration']
                        if 'Duration' in growthrun_sheet.columns
                        and not pd.isna(step['Duration'])
                        else 0 * ureg('minute').to('second').magnitude
                    ),
                    comment=f"Description: {typed_df_value(substrates_sheet, 'Description', str, substrate_index)}, Notes: {typed_df_value(substrates_sheet, 'Notes', str, substrate_index)}"
                    if 'Notes' and 'Description' in pregrowth_sheet.columns
                    else None,
                    environment=ChamberEnvironmentMovpe(
                        pressure=Pressure(
                            set_time=(pd.Series([0])),
                            set_value=pd.Series(
                                [
                                    (
                                        step['Chamber Pressure']
                                        if 'Chamber Pressure' in pregrowth_sheet.columns
                                        else None
                                    )
                                ]
                            )
                            * ureg('mbar').to('pascal').magnitude,
                        ),
                        rotation=Rotation(
                            set_time=(pd.Series([0])),
                            set_value=pd.Series(
                                [
                                    (
                                        step['Carrier Rotation']
                                        if 'Carrier Rotation' in pregrowth_sheet.columns
                                        else None
                                    )
                                ]
                            )
                            * ureg('rpm').to('rpm').magnitude,
                        ),
                        gas_flow=[
                            GasFlowMovpe(
                                gas=PubChemPureSubstanceSection(
                                    name=(
                                        step['Carrier Gas']
                                        if 'Carrier Gas' in pregrowth_sheet.columns
                                        else None
                                    ),
                                ),
                            ),
                        ],
                        uniform_gas_flow_rate=VolumetricFlowRate(
                            set_time=(pd.Series([0])),
                            set_value=pd.Series(
                                [
                                    (
                                        step['Carrier Gas Flow']
                                        if 'Carrier Gas Flow' in pregrowth_sheet.columns
                                        else None
                                    )
                                ]
                            )
                            * ureg('cm ** 3 / minute')
                            .to('meter ** 3 / second')
                            .magnitude,
                        ),
                    ),
                    sample_parameters=[
                        SampleParametersMovpe(
                            filament_temperature=FilamentTemperature(
                                time=(pd.Series([0])),
                                value=(
                                    pd.Series(
                                        [
                                            step['Substrate Temperature']
                                            if 'Substrate Temperature'
                                            in pregrowth_sheet.columns
                                            else None
                                        ]
                                    )
                                    * ureg('celsius').to('kelvin').magnitude
                                ),
                                set_time=(pd.Series([0])),
                                set_value=(
                                    pd.Series(
                                        [
                                            step['Substrate Temperature']
                                            if 'Substrate Temperature'
                                            in pregrowth_sheet.columns
                                            else None
                                        ]
                                    )
                                    * ureg('celsius').to('kelvin').magnitude
                                ),
                            )
                        )
                    ],
                )
            )

        # creating PRE-growth process objects
        pregrowth_process_object = GrowthMovpeIMEM(
            name='Pregrowth MOVPE',
            lab_id=sample_id,
            susceptor=typed_df_value(
                substrates_sheet, 'Susceptor', str, substrate_index
            ),
            mask=typed_df_value(substrates_sheet, 'Mask', str, substrate_index),
            pocket=typed_df_value(substrates_sheet, 'Pocket', str, substrate_index),
        )
        pregrowth_process_object.steps = pre_process_steps_lists

        # creating pregrowth process archives
        pregrowth_process_filename = (
            f'{sample_id}-pregrowth.GrowthMovpeIMEM.archive.{filetype}'
        )
        # Activity.normalize(pregrowth_process_object, archive, logger)
        pregrowth_process_archive = EntryArchive(
            data=pregrowth_process_object,
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
            afm_data = AFMmeasurement(
                name=str(
                    row['Sample']
                    if 'Sample' in characterization_sheet.columns and row['Sample']
                    else None
                )
                + f' afm {index}',
                datetime=(
                    row['Date']
                    if 'Date' in characterization_sheet.columns and row['Date']
                    else None
                ),
                samples=[
                    CompositeSystemReference(
                        lab_id=str(
                            row['Sample']
                            if 'Sample' in hrxrd_sheet.columns and row['Sample']
                            else characterization_sheet
                        ),
                    )
                ],
                results=[
                    AFMresults(
                        name=(
                            row['Notes']
                            if 'Notes' in characterization_sheet.columns
                            and row['Notes']
                            else None
                        ),
                        roughness=(
                            row['Roughness']
                            if 'Roughness' in characterization_sheet.columns
                            and row['Roughness']
                            else None
                        ),
                        surface_features=(
                            row['Surface Features']
                            if 'Surface Features' in characterization_sheet.columns
                            and row['Surface Features']
                            else None
                        ),
                    )
                ],
            )
            afm_filename = f'{row["Sample"]}.AFMmeasurement.archive.{filetype}'
            afm_archive = EntryArchive(
                data=afm_data,
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
                    name=str(
                        row['Sample']
                        if 'Sample' in hrxrd_sheet.columns and row['Sample']
                        else None
                    )
                    + f' afm {index}',
                    reference=f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.m_context.upload_id, afm_filename)}#data',
                )
            )

        # creating XRD object(s)
        xrd_measurements = []
        for index, row in hrxrd_sheet.iterrows():
            xrd_measurements.append(
                XRDmeasurementReference(
                    name=str(
                        row['Sample']
                        if 'Sample' in hrxrd_sheet.columns and row['Sample']
                        else None
                    )
                    + f' afm {index}',
                    sample_id=(
                        row['Sample']
                        if 'Sample' in hrxrd_sheet.columns and row['Sample']
                        else None
                    ),
                    phase=(
                        row['Phase']
                        if 'Phase' in hrxrd_sheet.columns and row['Phase']
                        else None
                    ),
                    peak_position_2theta=(
                        row['Peak Position - 2theta']
                        if 'Peak Position - 2theta' in hrxrd_sheet.columns
                        and row['Peak Position - 2theta']
                        else None
                    ),
                    peak_fwhm_2theta=(
                        row['Peak FWHM - 2theta']
                        if 'Peak FWHM - 2theta' in hrxrd_sheet.columns
                        and row['Peak FWHM - 2theta']
                        else None
                    ),
                    peak_position_omega=(
                        row['Peak Position - Omega']
                        if 'Peak Position - Omega' in hrxrd_sheet.columns
                        and row['Peak Position - Omega']
                        else None
                    ),
                    peak_fwhm_rocking_curve=(
                        row['Peak FWHM Rocking Curve']
                        if 'Peak FWHM Rocking Curve' in hrxrd_sheet.columns
                        and row['Peak FWHM Rocking Curve']
                        else None
                    ),
                    reflection=(
                        row['Reflection']
                        if 'Reflection' in hrxrd_sheet.columns and row['Reflection']
                        else None
                    ),
                    description=(
                        row['Notes']
                        if 'Notes' in hrxrd_sheet.columns and row['Notes']
                        else None
                    ),
                )
            )

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
            data=experiment_data,
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
