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
import re
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
    CompositeSystemReference,
    EntityReference,
    ElementalComposition,
    PureSubstanceComponent,
    PureSubstanceSection,
    PureSubstance,
    PureSubstanceSection,
    ExperimentStep,
)

from nomad_material_processing.general import (
    SubstrateReference,
    ThinFilmReference,
    ThinFilmStackReference,
    Parallelepiped,
    SubstrateCrystalProperties,
    Miscut,
    CartesianMiscut,
    ProjectedMiscutOrientation,
    Dopant,
    MillerIndices,
    CrystallographicDirection,
)
from nomad_material_processing.vapor_deposition.general import (
    Pressure,
    VolumetricFlowRate,
    GasFlow,
    Temperature,
)

from nomad_material_processing.vapor_deposition.cvd.general import (
    PartialVaporPressure,
    BubblerEvaporator,
    Rotation,
    BubblerSource,
    GasCylinderSource,
    GasCylinderEvaporator,
    PushPurgeGasFlow,
    MistSource,
    MistEvaporator,
    ComponentConcentration,
)

from imem_nomad_plugin.general.schema import (
    SampleCutIMEM,
)

from imem_nomad_plugin.characterization.schema import (
    AFMmeasurement,
    AFMresults,
    HallMeasurement,
    HallResults,
    ReflectanceMeasurement,
    ReflectanceResults,
    SEMmeasurement,
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
    ChamberEnvironmentMovpe,
    SubstrateMovpe,
    Shape,
    PureSubstanceIMEM,
    PrecursorReference,
    PressureIMEM,
    TemperatureIMEM,
    SampleParametersMovpe,
    FilamentTemperature,
    XRDmeasurementReference,
    AFMmeasurementReference,
    HallMeasurementReference,
    ReflectanceReference,
    SEMReference,
    UVAbsorbanceReference,
    CharacterizationMovpeIMEM,
    ContactSputtering,
    ContactSputteringReference,
)

from imem_nomad_plugin.utils import (
    create_archive,
    # fetch_substrate,
    import_excel_sheet,
    fill_quantity,
    rename_block_cols,
    clean_col_names,
    split_list_by_element,
)


def populate_bubblers(growth_step: pd.DataFrame, bub_first_col: str):
    """
    Populate dopants from the substrate row.

    The function usage implies that the repeated columns are named
    through the rename_block_cols function as follows:
    'Bubbler Material.1', 'Bubbler Temp.1', 'Bubbler Material.2', 'Bubbler Temp.2', ...
    and not as pandas default:
    'Bubbler Material', 'Bubbler Temp', 'Bubbler Material.1', 'Bubbler Temp.1', ...
    """
    groups = split_list_by_element(growth_step.index, bub_first_col)
    bubblers = []
    for columns in groups:
        if any(bub_first_col in item for item in columns):
            quantity, sep, index = columns[0].rpartition('.')
            bubbler_obj = BubblerSource()
            bubbler_obj.material = [PureSubstanceComponent()]
            bubbler_obj.material[0].pure_substance = PureSubstanceSection()
            bubbler_obj.vapor_source = BubblerEvaporator()
            bubbler_obj.vapor_source.temperature = TemperatureIMEM()
            bubbler_obj.vapor_source.pressure = PressureIMEM()
            bubbler_obj.vapor_source.total_flow_rate = VolumetricFlowRate()

            material = fill_quantity(growth_step, f'Bubbler Material{sep}{index}')
            bubbler_obj.name = material
            bubbler_obj.material[0].substance_name = material
            bubbler_obj.material[0].pure_substance.name = material

            valve = fill_quantity(growth_step, f'Bubbler Valve{sep}{index}')
            if valve:
                bubbler_obj.valve = valve

            bub_temp = fill_quantity(
                growth_step,
                f'Bubbler Temp{sep}{index}',
                read_unit='celsius',
                array=True,
            )
            if bub_temp is not None:
                bubbler_obj.vapor_source.temperature.set_time = [0]
                bubbler_obj.vapor_source.temperature.set_value = bub_temp

            bub_press = fill_quantity(
                growth_step,
                f'Bubbler Pressure{sep}{index}',
                read_unit='mbar',
                array=True,
            )
            if bub_press:
                bubbler_obj.vapor_source.pressure.set_time = [0]
                bubbler_obj.vapor_source.pressure.set_value = bub_press

            bub_flow = fill_quantity(
                growth_step,
                f'Bubbler Flow{sep}{index}',
                read_unit='cm ** 3 / minute',
                array=True,
            )
            if bub_flow:
                bubbler_obj.vapor_source.total_flow_rate.set_time = [0]
                bubbler_obj.vapor_source.total_flow_rate.set_value = bub_flow

            bub_dilution = fill_quantity(growth_step, f'Bubbler Dilution{sep}{index}')
            if bub_dilution:
                bubbler_obj.vapor_source.dilution = bub_dilution

            bubblers.append(bubbler_obj)
    return bubblers
    # precursor_partial_pressure=PartialVaporPressure(
    #     set_value=pd.Series(
    #         [
    #             growth_run_file.get(
    #                 f"Partial Pressure{'' if i == 0 else '.' + str(i)}",
    #                 0,
    #             )[line_number]
    #         ]
    #     ),
    # ),
    # vapor_molar_flow_rate=MolarFlowRate(
    #     set_value=pd.Series(
    #         [
    #             growth_run_file.get(
    #                 f"Bubbler Molar Flux{'' if i == 0 else '.' + str(i)}",
    #                 0,
    #             )[line_number]
    #         ]
    #     )
    #     * ureg('mol / minute').to('mol / second').magnitude,
    # ),


def populate_cylinder(growth_step: pd.DataFrame, gas_first_col: str):
    """
    Populate dopants from the substrate row.

    The function usage implies that the repeated columns are named
    through the rename_block_cols function as follows:
    'Gas Cylinder Material.1', 'Cylinder Pressure.1', 'Gas Cylinder Material.2', 'Cylinder Pressure.2', ...
    and not as pandas default:
    'Gas Cylinder Material', 'Cylinder Pressure', 'Gas Cylinder Material.1', 'Cylinder Pressure.1', ...
    """
    groups = split_list_by_element(growth_step.index, gas_first_col)
    cylinders = []
    for columns in groups:
        if any(gas_first_col in item for item in columns):
            quantity, sep, index = columns[0].rpartition('.')
            cylinder_obj = GasCylinderSource()
            cylinder_obj.material = [PureSubstanceComponent()]
            cylinder_obj.material[0].pure_substance = PureSubstanceSection()
            cylinder_obj.vapor_source = GasCylinderEvaporator()
            cylinder_obj.vapor_source.pressure = PressureIMEM()
            cylinder_obj.vapor_source.total_flow_rate = VolumetricFlowRate()
            cylinder_obj.vapor_source.effective_flow_rate = VolumetricFlowRate()

            material = fill_quantity(growth_step, f'Gas Cylinder Material{sep}{index}')
            cylinder_obj.name = material
            cylinder_obj.material[0].substance_name = material
            cylinder_obj.material[0].pure_substance.name = material

            valve = fill_quantity(growth_step, f'Gas Valve{sep}{index}')
            if valve:
                cylinder_obj.valve = valve

            dilution = fill_quantity(
                growth_step, f'Cylinder Pressure{sep}{index}', read_unit='mbar'
            )
            if dilution:
                cylinder_obj.vapor_source.dilution_in_cylinder = dilution

            pressure = fill_quantity(
                growth_step,
                f'Cylinder Pressure{sep}{index}',
                read_unit='mbar',
                array=True,
            )
            if pressure:
                cylinder_obj.vapor_source.pressure.set_time = [0]
                cylinder_obj.vapor_source.pressure.set_value = pressure

            flow = fill_quantity(
                growth_step,
                f'Flow from MFC{sep}{index}',
                read_unit='cm **3 / minute',
                array=True,
            )
            if flow:
                cylinder_obj.vapor_source.total_flow_rate.set_time = [0]
                cylinder_obj.vapor_source.total_flow_rate.set_value = flow

            effective_flow = fill_quantity(
                growth_step,
                f'Effective  Flow{sep}{index}',
                read_unit='cm **3 / minute',
                array=True,
            )
            if effective_flow:
                cylinder_obj.vapor_source.total_flow_rate.set_time = [0]
                cylinder_obj.vapor_source.total_flow_rate.set_value = effective_flow
            cylinders.append(cylinder_obj)
    return cylinders


def populate_mist_source(
    growth_step: pd.DataFrame,
    growthstep_first_col: str,
    mist_frame: pd.DataFrame,
    mist_first_col: str,
):
    """
    Populate dopants from the substrate row.

    The function usage implies that the repeated columns are named
    through the rename_block_cols function as follows:
    'Mist Material.1', 'Mist Pressure.1', 'Mist Material.2', 'Mist Pressure.2', ...
    and not as pandas default:
    'Mist Material', 'Mist Pressure', 'Mist Material.1', 'Mist Pressure.1', ...
    """
    source_groups = split_list_by_element(growth_step.index, growthstep_first_col)
    material_groups = []
    if (
        mist_frame is not None
        and isinstance(mist_frame, pd.DataFrame)
        and not mist_frame.empty
    ):
        material_groups = split_list_by_element(mist_frame.index, mist_first_col)
    mists = []
    for source_col in source_groups:
        if any(growthstep_first_col in item for item in source_col):
            quantity, sep, index = source_col[0].rpartition('.')
            mist_obj = MistSource()
            mist_obj.material = []
            mist_obj.vapor_source = MistEvaporator()
            mist_obj.vapor_source.temperature = TemperatureIMEM()
            mist_obj.vapor_source.total_flow_rate = VolumetricFlowRate()

            if (
                mist_frame is not None
                and isinstance(mist_frame, pd.DataFrame)
                and not mist_frame.empty
            ):
                item = fill_quantity(mist_frame, 'Item')
                if item:
                    mist_obj.item = item

                stirring = fill_quantity(mist_frame, 'Time')
                if stirring:
                    mist_obj.stirring_time = stirring

                notes = fill_quantity(mist_frame, 'Notes')
                if notes:
                    mist_obj.description = notes

                temperature = fill_quantity(
                    mist_frame, 'Temperature', read_unit='celsius', array=True
                )
                if temperature is not None:
                    mist_obj.vapor_source.temperature.set_time = [0]
                    mist_obj.vapor_source.temperature.set_value = temperature

            name = fill_quantity(growth_step, f'MIST Source 1{sep}{index}')
            if name:
                mist_obj.name = name

            valve = fill_quantity(growth_step, f'MIST Valve{sep}{index}')
            if valve:
                mist_obj.valve = valve

            flow = fill_quantity(
                growth_step,
                f'MIST Flow MFC{sep}{index}',
                read_unit='cm **3 / minute',
                array=True,
            )
            if flow:
                mist_obj.vapor_source.total_flow_rate.set_time = [0]
                mist_obj.vapor_source.total_flow_rate.set_value = flow

            for columns in material_groups:
                if any(mist_first_col in item for item in columns):
                    quantity2, sep2, index2 = columns[0].rpartition('.')

                    material = fill_quantity(mist_frame, f'Material{sep2}{index2}')
                    conc = fill_quantity(mist_frame, f'Concentration{sep2}{index2}')
                    mist_obj.material.append(
                        ComponentConcentration(
                            substance_name=material,
                            pure_substance=PureSubstanceSection(
                                name=material,
                            ),
                            theoretical_concentration=conc,
                            effective_concentration=conc,
                        )
                    )
            mists.append(mist_obj)
    return mists


def populate_elements(substrate_row: pd.DataFrame):
    column = 'Elements'
    elements = []
    for columns in substrate_row.index:
        if column in columns:
            quantity, sep, index = column.rpartition('.')
            element_obj = ElementalComposition()

            el = (
                fill_quantity(substrate_row, f'{columns}{sep}{index}')
                if f'{columns}{sep}{index}' in columns
                else fill_quantity(substrate_row, f'{columns}')
            )
            if el:
                element_obj.element = el
            elements.append(element_obj)
    return elements


def populate_dopants(substrate_row: pd.DataFrame, first_col: str):
    """
    Populate dopants from the substrate row.

    The function usage implies that the repeated columns are named
    through the rename_block_cols function as follows:
    'Doping species.1', 'Doping Level.1', 'Doping species.2', 'Doping Level.2', ...
    and not as pandas default:
    'Doping species', 'Doping Level', 'Doping species.1', 'Doping Level.1', ...

    caveats:
    - use a loop to iterate over rows (.iterrows()) in the dataframe and call this function on each row

    """
    groups = split_list_by_element(substrate_row.index, first_col)

    dopants = []
    for columns in groups:
        if any(first_col in item for item in columns):
            quantity, sep, index = columns[0].rpartition('.')
            dopant_obj = Dopant()

            species = fill_quantity(substrate_row, f'Doping species{sep}{index}')
            if species:
                dopant_obj.element = species

            dop_level = fill_quantity(substrate_row, f'Doping Level{sep}{index}')
            if dop_level:
                dopant_obj.doping_level = dop_level
            dopants.append(dopant_obj)
    return dopants


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

        ################################ excel reading ################################

        # creates a new DataFrame with the stripped column names.
        # sheet = pd.read_excel(xlsx, 'Overview', comment='#', converters={'Sample': str})
        # overview_sheet = sheet.rename(columns=lambda x: x.strip())
        if 'Overview' in xlsx.sheet_names:
            overview_sheet = import_excel_sheet(
                xlsx, 'Overview', converters={'Sample': str}
            )
        if 'Substrate' in xlsx.sheet_names:
            substrates_sheet = import_excel_sheet(
                xlsx,
                'Substrate',
                converters={'Orientation': str, 'Off-cut Orientation': str},
            )
            substrate_cols = clean_col_names(substrates_sheet)
            dopant_quantities = [
                'Doping species',
                'Doping Level',
            ]
            dopant_first_col = 'Doping Level'
            substrate_cols = rename_block_cols(
                substrate_cols, dopant_quantities, dopant_first_col
            )
            element_quantities = [
                'Elements',
            ]
            element_first_col = element_quantities[0]
            substrate_cols = rename_block_cols(
                substrate_cols, element_quantities, element_first_col
            )
            substrates_sheet.columns = substrate_cols
        if 'GrowthRun' in xlsx.sheet_names:
            growthrun_sheet = import_excel_sheet(xlsx, 'GrowthRun')
            growthrun_cols = clean_col_names(growthrun_sheet)
            bubbler_cols = [
                'Bubbler Temp',
                'Bubbler Pressure',
                'Partial Pressure',
                'Bubbler Dilution',
                'Inject',
                'Bubbler Flow',
                'Bubbler Material',
                'Bubbler Valve',
            ]
            bub_first_col = 'Bubbler Material'
            growthrun_cols = rename_block_cols(
                growthrun_cols, bubbler_cols, bub_first_col
            )
            gas_source_quantities = [
                'Gas Cylinder Material',
                'Dilution in Cylinder',
                'Flow from MFC',
                'Effective  Flow',
                'Gas Partial Pressure',
                'Cylinder Pressure',
                'Gas Valve',
            ]
            gas_first_col = 'Gas Cylinder Material'
            growthrun_cols = rename_block_cols(
                growthrun_cols, gas_source_quantities, gas_first_col
            )
            mist_quantities = [
                'MIST Source 1',
                'MIST Flow MFC',
                'MIST Valve',
            ]
            mist_first_col1 = 'MIST Source 1'
            growthrun_cols = rename_block_cols(
                growthrun_cols, mist_quantities, mist_first_col1
            )
            growthrun_sheet.columns = growthrun_cols
        if 'Precursors' in xlsx.sheet_names:
            precursors_sheet = import_excel_sheet(
                xlsx, 'Precursors', converters={'CAS': str}
            )
        if 'Mist' in xlsx.sheet_names:
            mist_sheet = import_excel_sheet(xlsx, 'Mist')
            mist_cols = clean_col_names(mist_sheet)
            mist_quantities = [
                'Material',
                'Concentration',
            ]
            mist_first_col2 = 'Material'
            mist_cols = rename_block_cols(mist_cols, mist_quantities, mist_first_col2)
            mist_sheet.columns = mist_cols
        if 'Pregrowth' in xlsx.sheet_names:
            pregrowth_sheet = import_excel_sheet(xlsx, 'Pregrowth')
        if 'SampleCut' in xlsx.sheet_names:
            samplecut_sheet = import_excel_sheet(xlsx, 'SampleCut')
        if 'HRXRD' in xlsx.sheet_names:
            hrxrd_sheet = import_excel_sheet(xlsx, 'HRXRD', converters={'Sample': str})
        if 'AFMReflectanceSEM' in xlsx.sheet_names:
            characterization_sheet = import_excel_sheet(
                xlsx, 'AFMReflectanceSEM', converters={'Sample': str}
            )
        if 'ElectroOptical' in xlsx.sheet_names:
            electro_optical_sheet = import_excel_sheet(xlsx, 'ElectroOptical')
        if 'Contacts' in xlsx.sheet_names:
            contacts_sheet = import_excel_sheet(xlsx, 'Contacts')

        try:
            sample_id = fill_quantity(overview_sheet.iloc[0], 'Sample')
        except IndexError:
            sample_id = None
        ############################## end excel reading ##############################

        # creating Substrate archives
        substrate_index = None
        for (
            substrate_index,
            substrate_row,
        ) in substrates_sheet.iterrows():
            substrate_id = fill_quantity(substrate_row, 'Substrates')
            # creating Substrate archives
            substrate_filename = (
                f'{substrate_id}_{substrate_index}.SubstrateIMEM.archive.{filetype}'
            )
            substrate_data = SubstrateMovpe()
            substrate_data.geometry = Shape()
            substrate_data.crystal_properties = SubstrateCrystalProperties()
            substrate_data.crystal_properties.surface_orientation = CrystallographicDirection()
            substrate_data.crystal_properties.surface_orientation.hkl_reciprocal = MillerIndices()
            substrate_data.crystal_properties.miscut = Miscut()
            substrate_data.crystal_properties.miscut.cartesian_miscut = (
                CartesianMiscut()
            )
            substrate_data.crystal_properties.miscut.cartesian_miscut.reference_orientation = ProjectedMiscutOrientation()
            substrate_data.crystal_properties.miscut.cartesian_miscut.reference_orientation.hkl_reciprocal = MillerIndices()

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
                substrate_data.crystal_properties.surface_orientation.hkl_reciprocal.h_index = int(orientation[0])
                substrate_data.crystal_properties.surface_orientation.hkl_reciprocal.k_index = int(orientation[1])
                substrate_data.crystal_properties.surface_orientation.hkl_reciprocal.l_index = int(orientation[2])

            miscut_angle = fill_quantity(substrate_row, 'Off-cut')
            if miscut_angle:
                substrate_data.crystal_properties.miscut.cartesian_miscut.reference_orientation.angle = miscut_angle

            mo = fill_quantity(substrate_row, 'Off-cut Orientation')
            if mo:
                substrate_data.crystal_properties.miscut.cartesian_miscut.reference_orientation.hkl_reciprocal.h_index = int(mo[0])
                substrate_data.crystal_properties.miscut.cartesian_miscut.reference_orientation.hkl_reciprocal.k_index = int(mo[1])
                substrate_data.crystal_properties.miscut.cartesian_miscut.reference_orientation.hkl_reciprocal.l_index = int(mo[2])

            substrate_data.elemental_composition = populate_elements(substrate_row)

            substrate_data.dopants = populate_dopants(substrate_row, dopant_first_col)

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
        # if a substrate will be located in a different upload, implement the following:
        # fetch_substrate(archive, sample_id, substrate_id, logger)
        # sleep(1.5)
        substrate_ref = None
        if substrate_id and sample_id:
            substrate_ref = f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.m_context.upload_id, substrate_filename)}#data'
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

        # creating Contact Sputtering archive
        samplecut_data = SampleCutIMEM(
            name=str(sample_id) + ' sample cut',
            parent_sample=ThinFilmStackMovpeReference(
                reference=f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.m_context.upload_id, grown_sample_filename)}#data'
            ),
        )
        sputterings = []
        for sput_id, row in contacts_sheet.iterrows():
            sputtering_data = ContactSputtering()

            sputtering_sample_id = (
                fill_quantity(row, 'Sample Piece')
                if fill_quantity(row, 'Sample Piece')
                else f'{sample_id} sputtering {sput_id} '
            )
            sputtering_data.samples = []
            if sputtering_sample_id:
                sputtering_data.name = str(sputtering_sample_id)
                sputtering_data.lab_id = str(sputtering_sample_id)
                sputtering_data.samples.append(
                    CompositeSystemReference(
                        lab_id=str(sputtering_sample_id),
                    )
                )
            contact_type = fill_quantity(row, 'Contact Type')
            if contact_type:
                sputtering_data.contact_type = contact_type
            geometry = fill_quantity(row, 'Geometry')
            if geometry:
                sputtering_data.geometry = geometry

            sputtering_filename = (
                f'{sputtering_sample_id}.Sputtering.archive.{filetype}'
            )
            sputtering_archive = EntryArchive(
                data=sputtering_data if sputtering_data else ContactSputtering(),
                m_context=archive.m_context,
                metadata=EntryMetadata(upload_id=archive.m_context.upload_id),
            )
            create_archive(
                sputtering_archive.m_to_dict(),
                archive.m_context,
                sputtering_filename,
                filetype,
                logger,
            )
            sputterings.append(
                ContactSputteringReference(
                    reference=f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.m_context.upload_id, sputtering_filename)}#data'
                )
            )

        # creating growth process step objects
        process_steps_lists = []
        for step_index, step in growthrun_sheet.iterrows():
            growth_step = GrowthStepMovpeIMEM()
            growth_step.environment = ChamberEnvironmentMovpe()
            growth_step.environment.pressure = PressureIMEM()
            growth_step.environment.rotation = Rotation()
            growth_step.environment.gas_flow = []
            growth_step.environment.uniform_gas = GasFlow()
            growth_step.environment.uniform_gas.flow_rate = VolumetricFlowRate()
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
                step, 'Pressure', read_unit='mbar', array=True
            )
            if growth_step_environment_pressure_set_value:
                growth_step.environment.pressure.set_time = [0]
                growth_step.environment.pressure.set_value = (
                    growth_step_environment_pressure_set_value
                )

            growth_step_environment_rotation_set_value = fill_quantity(
                step, 'Rotation', read_unit='rpm', array=True
            )
            if growth_step_environment_rotation_set_value:
                growth_step.environment.rotation.set_time = [0]
                growth_step.environment.rotation.set_value = (
                    growth_step_environment_rotation_set_value
                )

            growth_step_uniform_gas_flow_gas_name = fill_quantity(
                step, 'Carrier Gas'
            )
            if growth_step_uniform_gas_flow_gas_name:
                growth_step.environment.uniform_gas.gas=PureSubstanceSection(
                            name=growth_step_uniform_gas_flow_gas_name,
                        )
            growth_step_environment_uniform_gas_flow_rate_set_value = fill_quantity(
                step, 'Uniform Valve', read_unit='cm ** 3 / minute', array=True
            )
            if growth_step_environment_uniform_gas_flow_rate_set_value:
                growth_step.environment.uniform_gas.flow_rate.set_time = [0]
                growth_step.environment.uniform_gas.flow_rate.set_value = (
                    growth_step_environment_uniform_gas_flow_rate_set_value
                )

            growth_step_sample_parameters_filament_temperature_set_value = (
                fill_quantity(step, 'Temperature', read_unit='celsius', array=True)
            )
            if growth_step_sample_parameters_filament_temperature_set_value is not None:
                growth_step.sample_parameters[0].filament_temperature.set_time = [0]
                growth_step.sample_parameters[
                    0
                ].filament_temperature.set_value = (
                    growth_step_sample_parameters_filament_temperature_set_value
                )

            # the assumption stands here that the is only ONE MIST SOURCE
            mist_row = None
            for mist_index, row in mist_sheet.iterrows():
                if mist_index == 0:  # <--------------------------
                    mist_row = row
                    break

            bubblers = populate_bubblers(step, bub_first_col)
            gas_cylinders = populate_cylinder(step, gas_first_col)
            mist = populate_mist_source(
                step, mist_first_col1, mist_row, mist_first_col2
            )

            sources = []
            if bubblers:
                sources += bubblers
            if gas_cylinders:
                sources += gas_cylinders
            if mist:
                sources += mist

            growth_step.sources = sources
            process_steps_lists.append(growth_step)

        # creating growth process objects
        growth_process_object = GrowthMovpeIMEM()

        growth_process_object.name = f'{sample_id} growth'
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

        # creating growth PRE-process step objects
        pre_process_steps_lists = []
        for step_id, step in pregrowth_sheet.iterrows():
            pregrowth_step = GrowthStepMovpeIMEM()
            pregrowth_step.environment = ChamberEnvironmentMovpe()
            pregrowth_step.environment.pressure = PressureIMEM()
            pregrowth_step.environment.rotation = Rotation()
            pregrowth_step.environment.gas_flow = []
            pregrowth_step.environment.uniform_gas = GasFlow()
            pregrowth_step.environment.uniform_gas.flow_rate = VolumetricFlowRate()
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
                step, 'Chamber Pressure', read_unit='mbar', array=True
            )
            if pregrowth_step_environment_pressure_set_value:
                pregrowth_step.environment.pressure.set_time = [0]
                pregrowth_step.environment.pressure.set_value = (
                    pregrowth_step_environment_pressure_set_value
                )

            pregrowth_step_environment_rotation_set_value = fill_quantity(
                step, 'Carrier Rotation', read_unit='rpm', array=True
            )
            if pregrowth_step_environment_rotation_set_value:
                pregrowth_step.environment.rotation.set_time = [0]
                pregrowth_step.environment.rotation.set_value = (
                    pregrowth_step_environment_rotation_set_value
                )

            pregrowth_step_uniform_gas_flow_gas_name = fill_quantity(
                step, 'Carrier Gas'
            )
            if pregrowth_step_uniform_gas_flow_gas_name:
                pregrowth_step.environment.uniform_gas.gas=PureSubstanceSection(
                            name=pregrowth_step_uniform_gas_flow_gas_name,
                        )
            pregrowth_step_uniform_gas_flow_rate_set_value = fill_quantity(
                step, 'Carrier Gas Flow', read_unit='cm ** 3 / minute', array=True
            )
            if pregrowth_step_uniform_gas_flow_rate_set_value:
                pregrowth_step.environment.uniform_gas.flow_rate.set_time = [0]
                pregrowth_step.environment.uniform_gas.flow_rate.set_value = (
                    pregrowth_step_uniform_gas_flow_rate_set_value
                )

            pregrowth_step_sample_parameters_filament_temperature_set_value = (
                fill_quantity(
                    step, 'Substrate Temperature', read_unit='celsius', array=True
                )
            )
            if (
                pregrowth_step_sample_parameters_filament_temperature_set_value
                is not None
            ):
                pregrowth_step.sample_parameters[0].filament_temperature.set_time = [0]
                pregrowth_step.sample_parameters[
                    0
                ].filament_temperature.set_value = (
                    pregrowth_step_sample_parameters_filament_temperature_set_value
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
        pregrowth_process_object.steps = pre_process_steps_lists

        # creating AFM archive
        afm_measurements = []
        for index, row in characterization_sheet.iterrows():
            afm_data = AFMmeasurement()
            afm_data.samples = []
            afm_data.results = [AFMresults()]

            afm_name = (
                fill_quantity(row, 'Sample')
                if fill_quantity(row, 'Sample')
                else sample_id
            )
            if afm_name:
                afm_data.name = f'{afm_name} afm {index}'
                afm_data.samples.append(
                    CompositeSystemReference(
                        lab_id=str(afm_name),
                    )
                )
                afm_filename = f'{afm_name}_{index}.AFM.archive.{filetype}'
            else:
                afm_data.name = f'{sample_id} afm {index}'
                afm_filename = f'{sample_id}_{index}.AFMm.archive.{filetype}'

            afm_datetime = fill_quantity(row, 'Date')
            if afm_datetime:
                afm_data.datetime = afm_datetime

            afm_results_notes = fill_quantity(row, 'Notes')
            if afm_results_notes:
                afm_data.results[0].name = afm_results_notes

            afm_results_roughness = fill_quantity(
                row, 'Roughness', read_unit='nanometer'
            )
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

        # creating Hall archive
        hall_measurements = []
        for index, row in electro_optical_sheet.iterrows():
            hall_data = HallMeasurement()
            hall_data.samples = []
            hall_data.results = [HallResults()]

            hall_name = (
                fill_quantity(row, 'Sample')
                if fill_quantity(row, 'Sample')
                else sample_id
            )
            if hall_name:
                hall_data.name = f'{hall_name} hall {index}'
                hall_data.samples.append(
                    CompositeSystemReference(
                        lab_id=str(hall_name),
                    )
                )
                hall_filename = f'{hall_name}_{index}.Hall.archive.{filetype}'
            else:
                hall_data.name = f'{sample_id} hall {index}'
                hall_filename = f'{sample_id}_{index}.Hall.archive.{filetype}'

            hall_datetime = fill_quantity(row, 'Date')
            if hall_datetime:
                hall_data.datetime = hall_datetime

            hall_results_notes = fill_quantity(row, 'Notes')
            if hall_results_notes:
                hall_data.results[0].name = hall_results_notes

            resistivity = fill_quantity(row, 'Resistivity', read_unit='ohm * cm')
            if resistivity:
                hall_data.results[0].resistivity = resistivity

            mobility = fill_quantity(row, 'Mobility', read_unit='cm ** 2 / V / s')
            if mobility:
                hall_data.results[0].hall_mobility = mobility

            carrier_concentration = fill_quantity(
                row, 'Carrier Concentration', read_unit='1 / cm ** 3'
            )
            if carrier_concentration:
                hall_data.results[0].carrier_density = carrier_concentration

            hall_archive = EntryArchive(
                data=hall_data if hall_data else HallMeasurement(),
                m_context=archive.m_context,
                metadata=EntryMetadata(upload_id=archive.m_context.upload_id),
            )
            create_archive(
                hall_archive.m_to_dict(),
                archive.m_context,
                hall_filename,
                filetype,
                logger,
            )
            hall_measurements.append(
                HallMeasurementReference(
                    name=f'{hall_name} hall {index}',
                    reference=f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.m_context.upload_id, hall_filename)}#data',
                )
            )

        # creating reflectance archive
        reflec_measurements = []
        for index, row in characterization_sheet.iterrows():
            reflec_data = ReflectanceMeasurement()
            reflec_data.samples = []
            reflec_data.results = [ReflectanceResults()]

            reflec_name = (
                fill_quantity(row, 'Sample')
                if fill_quantity(row, 'Sample')
                else sample_id
            )
            if reflec_name:
                reflec_data.name = f'{reflec_name} reflectance {index}'
                reflec_data.samples.append(
                    CompositeSystemReference(
                        lab_id=str(reflec_name),
                    )
                )
                reflec_filename = (
                    f'{reflec_name}_{index}.Reflectance.archive.{filetype}'
                )
            else:
                reflec_data.name = f'{sample_id} reflectance {index}'
                reflec_filename = f'{sample_id}_{index}.Reflectance.archive.{filetype}'

            datetime = fill_quantity(row, 'Date')
            if datetime:
                reflec_data.datetime = datetime

            results_notes = fill_quantity(row, 'Notes')
            if results_notes:
                reflec_data.results[0].name = results_notes

            thickness = fill_quantity(row, 'Thickness', read_unit='nanometer')
            if thickness:
                reflec_data.results[0].thickness = thickness

            x_coord = fill_quantity(row, 'X value')
            if x_coord:
                reflec_data.results[0].x_coordinate = x_coord

            y_coord = fill_quantity(row, 'Y value')
            if y_coord:
                reflec_data.results[0].y_coordinate = y_coord

            growth_rate = fill_quantity(
                row, 'Growth Rate', read_unit='nanometer/minute'
            )
            if growth_rate:
                reflec_data.results[0].growth_rate = growth_rate

            reflec_archive = EntryArchive(
                data=reflec_data if reflec_data else ReflectanceMeasurement(),
                m_context=archive.m_context,
                metadata=EntryMetadata(upload_id=archive.m_context.upload_id),
            )
            create_archive(
                reflec_archive.m_to_dict(),
                archive.m_context,
                reflec_filename,
                filetype,
                logger,
            )
            reflec_measurements.append(
                ReflectanceReference(
                    name=f'{reflec_name} reflectance {index}',
                    reference=f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.m_context.upload_id, reflec_filename)}#data',
                )
            )

        # creating SEM object(s)
        sem_measurements = []
        for index, row in characterization_sheet.iterrows():
            sem_measurement = SEMmeasurement()

            sem_name = (
                fill_quantity(row, 'Sample')
                if fill_quantity(row, 'Sample')
                else sample_id
            )
            if sem_name:
                sem_measurement.name = f'{sem_name} sem {index}'
                sem_filename = f'{reflec_name}_{index}.SEM.archive.{filetype}'
            else:
                sem_measurement.name = f'{sample_id} SEM {index}'
                sem_filename = f'{sample_id}_{index}.SEM.archive.{filetype}'

            sem_date = fill_quantity(row, 'Date')
            if sem_date:
                sem_measurement.datetime = sem_date

            surface = fill_quantity(row, 'Surface State')
            if surface:
                sem_measurement.surface_state = surface

            sem_archive = EntryArchive(
                data=sem_measurement if sem_measurement else SEMmeasurement(),
                m_context=archive.m_context,
                metadata=EntryMetadata(upload_id=archive.m_context.upload_id),
            )
            create_archive(
                sem_archive.m_to_dict(),
                archive.m_context,
                sem_filename,
                filetype,
                logger,
            )
            sem_measurements.append(
                SEMReference(
                    name=f'{sem_name} SEM {index}',
                    reference=f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.m_context.upload_id, sem_filename)}#data',
                )
            )

        # creating Absorbance object(s)
        abs_measurements = []
        for index, row in characterization_sheet.iterrows():
            abs_measurement = UVAbsorbanceReference()

            abs_name = (
                fill_quantity(row, 'Sample')
                if fill_quantity(row, 'Sample')
                else sample_id
            )
            if abs_name:
                abs_measurement.name = f'{abs_name} sem {index}'

            abs_date = fill_quantity(row, 'Date')
            if abs_date:
                abs_measurement.datetime = abs_date

            energy_gap = fill_quantity(row, 'Energy Gap ABS', read_unit='eV')
            if energy_gap:
                abs_measurement.energy_gap = energy_gap

            abs_coefficient = fill_quantity(row, 'ABS coefficient', read_unit='1/cm')
            if abs_coefficient:
                abs_measurement.abs_coefficient = abs_coefficient

            thickness = fill_quantity(row, 'Thickness', read_unit='nanometer')
            if thickness:
                abs_measurement.thickness = thickness

            abs_measurements.append(abs_measurement)

        # creating XRD object(s)
        xrd_measurements = []
        for index, row in hrxrd_sheet.iterrows():
            xrd_measurement = XRDmeasurementReference()

            xrd_name = (
                fill_quantity(row, 'Sample')
                if fill_quantity(row, 'Sample')
                else sample_id
            )
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

        # creating precursors archives
        precursors = []
        for _, row in precursors_sheet.iterrows():
            precursor = PureSubstanceIMEM()
            precursor.pure_substance = PureSubstanceSection()

            cas = fill_quantity(row, 'CAS')
            if cas:
                precursor.pure_substance.cas_number = str(cas)
            name = fill_quantity(row, 'Name')
            if name:
                precursor.name = name
                precursor_filename = f'{name}.PureSubstance.archive.{filetype}'
            description = fill_quantity(row, 'Description')
            notes = fill_quantity(row, 'Notes')
            if description:
                precursor.description = f'{description} - {notes}'
            buying_date = fill_quantity(row, 'Buying Date')
            if buying_date:
                precursor.datetime = buying_date
            supplier = fill_quantity(row, 'Supplier')
            if supplier:
                precursor.supplier = supplier
            opening_date = fill_quantity(row, 'Opening Date')
            if opening_date:
                precursor.opening_date = opening_date
            phase = fill_quantity(row, 'Phase')
            if phase:
                precursor.phase = phase
            pre_archive = EntryArchive(
                data=precursor if precursor else PureSubstanceIMEM(),
                m_context=archive.m_context,
                metadata=EntryMetadata(upload_id=archive.m_context.upload_id),
            )
            create_archive(
                pre_archive.m_to_dict(),
                archive.m_context,
                precursor_filename,
                filetype,
                logger,
            )
            precursors.append(
                PrecursorReference(
                    name=f'{name}',
                    reference=f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.m_context.upload_id, precursor_filename)}#data',
                )
            )
        
        growth_process_filename = (
            f'{sample_id}-growth.GrowthMovpeIMEM.archive.{filetype}'
        )
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
        )        # creating pregrowth process archives
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

        # EXPERIMENT
        # creating experiment archive
        sub_ref = None
        if substrate_ref is not None:
            sub_ref = SubstrateReference(reference=substrate_ref)
        else:
            sub_ref = SubstrateReference(name=substrate_id, lab_id=substrate_id)

        experiment_data = ExperimentMovpeIMEM(
            name=f'experiment {sample_id}',
            pregrowth=GrowthMovpeIMEMReference(
                reference=f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.m_context.upload_id, pregrowth_process_filename)}#data',
            ),
            growth_run=GrowthMovpeIMEMReference(
                reference=f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.m_context.upload_id, growth_process_filename)}#data',
            ),
            substrates=[sub_ref],
            sample_cut=SampleCutIMEMReference(
                reference=f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.m_context.upload_id, samplecut_filename)}#data',
            ),
            contact_sputtering=sputterings,
            characterization=CharacterizationMovpeIMEM(
                xrd=xrd_measurements,
                afm=afm_measurements,
                hall=hall_measurements,
                reflectance=reflec_measurements,
                sem=sem_measurements,
            ),
            precursors=precursors,
        )

        # add precursors to growth and pregrowth
        growth_process_object.precursors = precursors
        pregrowth_process_object.precursors = precursors

        # creating experiment archive

        experiment_filename = f'{sample_id}.ExperimentMovpeIMEM.archive.{filetype}'
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
