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

        sample_id = fill_quantity(overview_sheet, 'Sample', 0)

        # creating Substrate archives
        for substrate_index, substrate_id in enumerate(substrates_sheet['Substrates']):
            # creating Substrate archives
            substrate_filename = (
                f'{substrate_id}_{substrate_index}.SubstrateIKZ.archive.{filetype}'
            )
            substrate_data = SubstrateMovpe(
                lab_id=substrate_id,
                supplier_id=fill_quantity(
                    substrates_sheet, 'Substrate ID', substrate_index
                ),
                supplier=fill_quantity(substrates_sheet, 'Supplier', substrate_index),
                name=fill_quantity(substrates_sheet, 'Material', substrate_index),
                description=f"Description: {fill_quantity(substrates_sheet, 'Description', substrate_index)}, Notes: {fill_quantity(substrates_sheet, 'Notes', substrate_index)}",
                geometry=Shape(
                    width=fill_quantity(
                        substrates_sheet,
                        'Size X',
                        substrate_index,
                        read_unit='millimeter',
                    ),
                    length=fill_quantity(
                        substrates_sheet,
                        'Size Y',
                        substrate_index,
                        read_unit='millimeter',
                    ),
                    diameter=fill_quantity(
                        substrates_sheet,
                        'Size Diameter',
                        substrate_index,
                        read_unit='millimeter',
                    ),
                ),
                crystal_properties=SubstrateCrystalPropertiesMovpe(
                    orientation=(
                        fill_quantity(substrates_sheet, 'Orientation', substrate_index)
                    ),
                    miscut=MiscutMovpe(
                        angle=fill_quantity(
                            substrates_sheet, 'Off-cut', substrate_index
                        ),
                        orientation=fill_quantity(
                            substrates_sheet,
                            'Off-cut Orientation',
                            substrate_index,
                        ),
                    ),
                ),
                elemental_composition=populate_element(
                    substrate_index, substrates_sheet
                ),
                dopants=populate_dopant(substrate_index, substrates_sheet),
                annealing=fill_quantity(substrates_sheet, 'Annealing', substrate_index),
                cleaning=fill_quantity(substrates_sheet, 'Cleaning', substrate_index),
                regrowth=fill_quantity(substrates_sheet, 'Regrowth', substrate_index),
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
        ### TODO try to get rid of the fetch_substarte as it slows down the processing
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
                    name=fill_quantity(step, 'Name'),
                    step_index=step_index + 1,
                    duration=fill_quantity(step, 'Duration', read_unit='minute'),
                    comment=fill_quantity(step, 'Notes'),
                    sources=populate_sources(step_index, growthrun_sheet)
                    + populate_gas_source(step_index, growthrun_sheet),
                    environment=ChamberEnvironmentMovpe(
                        pressure=Pressure(
                            set_time=pd.Series([0]),
                            set_value=pd.Series(
                                [fill_quantity(step, 'Pressure', read_unit='mbar')]
                            ),
                        ),
                        rotation=Rotation(
                            set_time=pd.Series([0]),
                            set_value=pd.Series(
                                [fill_quantity(step, 'Rotation', read_unit='rpm')]
                            ),
                        ),
                        gas_flow=[
                            GasFlowMovpe(
                                gas=PubChemPureSubstanceSection(
                                    name=fill_quantity(step, 'Carrier Gas'),
                                ),
                            ),
                        ],
                        uniform_gas_flow_rate=VolumetricFlowRate(
                            set_time=pd.Series([0]),
                            set_value=pd.Series(
                                [
                                    fill_quantity(
                                        step,
                                        'Uniform Valve',
                                        read_unit='cm ** 3 / minute',
                                    )
                                ]
                            ),
                        ),
                    ),
                    sample_parameters=[
                        SampleParametersMovpe(
                            filament_temperature=FilamentTemperature(
                                time=pd.Series([0]),
                                value=pd.Series(
                                    [
                                        fill_quantity(
                                            step, 'Temperature', read_unit='celsius'
                                        )
                                    ]
                                ),
                                set_time=pd.Series([0]),
                                set_value=pd.Series(
                                    [
                                        fill_quantity(
                                            step, 'Temperature', read_unit='celsius'
                                        )
                                    ]
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
                susceptor=fill_quantity(substrates_sheet, 'Susceptor', substrate_index)
                if fill_quantity(substrates_sheet, 'Susceptor', substrate_index)
                else '-',
                mask=fill_quantity(substrates_sheet, 'Mask', substrate_index)
                if fill_quantity(substrates_sheet, 'Mask', substrate_index)
                else '-',
                pocket=fill_quantity(substrates_sheet, 'Pocket', substrate_index)
                if fill_quantity(substrates_sheet, 'Pocket', substrate_index)
                else '-',
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
                    name=fill_quantity(step, 'Step Name'),
                    step_index=step_id + 1,
                    duration=fill_quantity(step, 'Duration', read_unit='minute'),
                    comment=f"Description: {fill_quantity(substrates_sheet, 'Description', substrate_index)}, Notes: {fill_quantity(substrates_sheet, 'Notes', substrate_index)}",
                    environment=ChamberEnvironmentMovpe(
                        # pressure=Pressure(
                        #     set_time=pd.Series([0]),
                        #     set_value=pd.Series(
                        #         [
                        #             fill_quantity(
                        #                 step, 'Chamber Pressure', read_unit='mbar'
                        #             )
                        #         ]
                        #     ),
                        # ),
                        rotation=Rotation(
                            set_time=pd.Series([0]),
                            set_value=pd.Series(
                                [
                                    fill_quantity(
                                        step, 'Carrier Rotation', read_unit='rpm'
                                    )
                                ]
                            ),
                        ),
                        gas_flow=[
                            GasFlowMovpe(
                                gas=PubChemPureSubstanceSection(
                                    name=fill_quantity(step, 'Carrier Gas'),
                                ),
                            ),
                        ],
                        uniform_gas_flow_rate=VolumetricFlowRate(
                            set_time=pd.Series([0]),
                            set_value=pd.Series(
                                [
                                    fill_quantity(
                                        step,
                                        'Carrier Gas Flow',
                                        read_unit='cm ** 3 / minute',
                                    )
                                ]
                            ),
                        ),
                    ),
                    sample_parameters=[
                        SampleParametersMovpe(
                            filament_temperature=FilamentTemperature(
                                time=pd.Series([0]),
                                value=pd.Series(
                                    [
                                        fill_quantity(
                                            step,
                                            'Substrate Temperature',
                                            read_unit='celsius',
                                        )
                                    ]
                                ),
                                set_time=pd.Series([0]),
                                set_value=pd.Series(
                                    [
                                        fill_quantity(
                                            step,
                                            'Substrate Temperature',
                                            read_unit='celsius',
                                        )
                                    ]
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
            susceptor=fill_quantity(substrates_sheet, 'Susceptor', substrate_index),
            mask=fill_quantity(substrates_sheet, 'Mask', substrate_index),
            pocket=fill_quantity(substrates_sheet, 'Pocket', substrate_index),
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
                name=str(fill_quantity(row, 'Sample')) + f' afm {index}',
                datetime=fill_quantity(row, 'Date'),
                samples=[
                    CompositeSystemReference(
                        lab_id=str(fill_quantity(row, 'Sample')),
                    )
                ],
                results=[
                    AFMresults(
                        name=fill_quantity(row, 'Notes'),
                        roughness=fill_quantity(row, 'Roughness'),
                        surface_features=fill_quantity(row, 'Surface Features'),
                    )
                ],
            )
            afm_filename = (
                f'{fill_quantity(row,"Sample")}.AFMmeasurement.archive.{filetype}'
            )
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
                    name=fill_quantity(row, 'Sample') + f' afm {index}',
                    reference=f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.m_context.upload_id, afm_filename)}#data',
                )
            )

        # creating XRD object(s)
        xrd_measurements = []
        for index, row in hrxrd_sheet.iterrows():
            xrd_measurements.append(
                XRDmeasurementReference(
                    name=fill_quantity(row, 'Sample') + f' xrd {index}',
                    sample_id=fill_quantity(row, 'Sample') + f' xrd {index}',
                    phase=fill_quantity(row, 'Phase') + f' xrd {index}',
                    peak_position_2theta=fill_quantity(row, 'Peak Position - 2theta'),
                    peak_fwhm_2theta=fill_quantity(row, 'Peak FWHM - 2theta'),
                    peak_position_omega=fill_quantity(row, 'Peak Position - Omega'),
                    peak_fwhm_rocking_curve=fill_quantity(
                        row, 'Peak FWHM Rocking Curve'
                    ),
                    reflection=fill_quantity(row, 'Reflection'),
                    description=fill_quantity(row, 'Notes'),
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
