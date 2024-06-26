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

import re
from time import sleep, perf_counter
import pandas as pd
from typing import Dict, List
import yaml
import json
import math

from nomad.datamodel.context import ClientContext
from nomad.datamodel import EntryArchive
from nomad.metainfo import MSection, Quantity, Section
from nomad.parsing import MatchingParser
from nomad.datamodel.metainfo.annotations import ELNAnnotation
from nomad.datamodel.data import EntryData
from nomad.units import ureg
from nomad.datamodel.metainfo.basesections import (
    SystemComponent,
    CompositeSystemReference,
    ElementalComposition,
    PureSubstanceComponent,
    PureSubstanceSection,
    ExperimentStep,
)

from nomad_material_processing import Dopant, SubstrateReference
from nomad_material_processing.vapor_deposition import (
    MolarFlowRate,
    Temperature,
    Pressure,
    VolumetricFlowRate,
)
from nomad_material_processing.vapor_deposition.cvd import (
    PartialVaporPressure,
    BubblerEvaporator,
)

from nomad.datamodel.datamodel import EntryArchive, EntryMetadata


def get_reference(upload_id, entry_id):
    return f'../uploads/{upload_id}/archive/{entry_id}'


def get_entry_id(upload_id, filename):
    from nomad.utils import hash

    return hash(upload_id, filename)


def get_hash_ref(upload_id, filename):
    return f'{get_reference(upload_id, get_entry_id(upload_id, filename))}#data'


def nan_equal(a, b):
    """
    Compare two values with NaN values.
    """
    if isinstance(a, float) and isinstance(b, float):
        return a == b or (math.isnan(a) and math.isnan(b))
    elif isinstance(a, dict) and isinstance(b, dict):
        return dict_nan_equal(a, b)
    elif isinstance(a, list) and isinstance(b, list):
        return list_nan_equal(a, b)
    else:
        return a == b


def list_nan_equal(list1, list2):
    """
    Compare two lists with NaN values.
    """
    if len(list1) != len(list2):
        return False
    for a, b in zip(list1, list2):
        if not nan_equal(a, b):
            return False
    return True


def dict_nan_equal(dict1, dict2):
    """
    Compare two dictionaries with NaN values.
    """
    if set(dict1.keys()) != set(dict2.keys()):
        return False
    for key in dict1:
        if not nan_equal(dict1[key], dict2[key]):
            return False
    return True


def create_archive(
    entry_dict, context, filename, file_type, logger, *, overwrite: bool = False
):
    from nomad.datamodel.context import ClientContext
    from nomad.datamodel import EntryArchive

    file_exists = context.raw_path_exists(filename)
    dicts_are_equal = None
    if isinstance(context, ClientContext):
        return None
    if file_exists:
        with context.raw_file(filename, 'r') as file:
            existing_dict = yaml.safe_load(file)
            dicts_are_equal = dict_nan_equal(existing_dict, entry_dict)
    if not file_exists or overwrite or dicts_are_equal:
        with context.raw_file(filename, 'w') as newfile:
            if file_type == 'json':
                json.dump(entry_dict, newfile)
            elif file_type == 'yaml':
                yaml.dump(entry_dict, newfile)
        context.upload.process_updated_raw_file(filename, allow_modify=True)
    elif file_exists and not overwrite and not dicts_are_equal:
        logger.error(
            f'{filename} archive file already exists. '
            f'You are trying to overwrite it with a different content. '
            f'To do so, remove the existing archive and click reprocess again.'
        )
    return get_hash_ref(context.upload_id, filename)

    # !! useful to fetch the upload_id from another upload.
    # experiment_context = ServerContext(
    #         get_upload_with_read_access(
    #             matches["upload_id"][0],
    #             User(
    #                 is_admin=True,
    #                 user_id=current_parse_archive.metadata.main_author.user_id,
    #             ),
    #             include_others=True,
    #         )
    #     )  # Upload(upload_id=matches["upload_id"][0]))


# def create_archive(
#     entry_dict, context, file_name, file_type, logger, *, bypass_check: bool = False
# ):
#     import yaml
#     import json

#     if not context.raw_path_exists(file_name) or bypass_check:
#         with context.raw_file(file_name, "w") as outfile:
#             if file_type == "json":
#                 json.dump(entry_dict, outfile)
#             elif file_type == "yaml":
#                 yaml.dump(entry_dict, outfile)
#         context.upload.process_updated_raw_file(file_name, allow_modify=True)
#     else:
#         logger.error(
#             f"{file_name} archive file already exists."
#             f"If you intend to reprocess the older archive file, remove the existing one and run reprocessing again."
#         )


def df_value(dataframe, column_header, index=None):
    """
    Fetches a value from a DataFrame.
    """
    if column_header in dataframe.columns:
        if index is not None:
            return dataframe[column_header][index]
        return dataframe[column_header]
    return None


def typed_df_value(dataframe, column_header, value_type, index=None):
    """
    Fetches a value of a specified type from a DataFrame.
    """
    value = df_value(dataframe, column_header, index)
    if value_type is str:
        return str(value)
    else:
        return value


def row_to_array(dataframe: pd.DataFrame, quantity: str, row_index: int) -> pd.Series:
    """
    Extracts values from a DataFrame row across multiple columns with similar names.

    This function takes a DataFrame, a list of header names, and a row index. It extracts
    the values from the specified row across all columns whose names start with any of the
    specified header names. The column names are expected to follow a specific pattern:
    the base header name, followed by a dot and an integer index (e.g., 'header.1', 'header.2', etc.).
    The function returns a Series containing all the extracted values.

    Args:
        dataframe (pd.DataFrame): The DataFrame containing the data.
        quantities (List[str]): The base names of the headers. The DataFrame should contain
                                 multiple columns with names that start with each base name.
        row_index (int): The index of the row from which to extract the values.

    Returns:
        pd.Series: A Series containing all the values extracted from the specified row across
                   all columns with names that start with any of the specified header names.

    Example:
        >>> df = pd.DataFrame({
        ...     'header': [1, 2, 3],
        ...     'header.1': [4, 5, 6],
        ...     'header.2': [7, 8, 9]
        ... })
        >>> array = row_to_array(df, ['header'], 0)
        >>> print(array)
        0    1
        1    4
        2    7
        dtype: int64
    """
    column_names = [col for col in dataframe.columns if col.startswith(quantity)]
    array = pd.Series(
        [
            dataframe[col].loc[row_index]
            for col in column_names
            if col is float or col is int
        ]
    )
    return array


def clean_timeseries(time_array: pd.Series, value_array: pd.Series) -> pd.Series:
    """
    clean time and value array pairs by removing NaNs
    """
    # for i in time_array.index:
    #     if pd.isna(time_array.loc[i]).any() or pd.isna(value_array.loc[i]).any():
    #         time_array = time_array.drop(i)
    #         value_array = value_array.drop(i)
    # return time_array, value_array
    df = pd.concat([time_array, value_array], axis=1)
    df = df.dropna()
    return df.iloc[:, 1], df.iloc[:, 0]


def row_timeseries(
    dataframe: pd.DataFrame, time_header: str, value_header: str, row_index: int
) -> pd.Series:
    """
    Extracts and cleans a timeseries from a row of a DataFrame.

    This function takes a DataFrame, two header names (one for time and one for values),
    and a row index. It extracts the timeseries data from the specified row, where the
    time and value data are stored in columns with the specified header names. The function
    then cleans the extracted timeseries by removing any NaN values.

    Args:
        dataframe (pd.DataFrame): The DataFrame containing the timeseries data.
        time_header (str): The header name for the time data. The DataFrame should contain
                            multiple columns with this header name, each containing a time point.
        value_header (str): The header name for the value data. The DataFrame should contain
                             multiple columns with this header name, each containing a value
                             corresponding to a time point.
        row_index (int): The index of the row from which to extract the timeseries.

    Returns:
        tuple: A tuple containing two pd.Series objects. The first Series contains the time
               points of the timeseries, and the second Series contains the corresponding values.
               Both Series have been cleaned to remove any NaN values.

    Example:
        >>> df = pd.DataFrame({
        ...     'time': [1, 2, 3],
        ...     'time.1': [4, 5, 6],
        ...     'value': [7, 8, 9],
        ...     'value.1': [10, 11, 12]
        ... })
        >>> time_series, value_series = row_timeseries(df, 'time', 'value', 0)
        >>> print(time_series)
        0    1
        1    4
        dtype: int64
        >>> print(value_series)
        0     7
        1    10
        dtype: int64
    """
    return clean_timeseries(
        row_to_array(
            dataframe,
            value_header,
            row_index,
        ),
        row_to_array(
            dataframe,
            time_header,
            row_index,
        ),
    )


def clean_dataframe_headers(dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Picks first row of a DataFrame and makes it as new headers.

    This function takes a DataFrame (without headers) as input and uses the first row as the new column names.
    It modifies these column names by removing trailing spaces and handling duplicate column names.
    Duplicate column names are handled by appending a count to the end of the column name.
    This functionality also a built-in in pandas, but it is not used here because it does not handle
    duplicate column names as expected.

    After setting the new column names, the function removes the first row (which has been used as headers)
    and resets the index.

    Args:
        dataframe (pd.DataFrame): The input DataFrame (without headers) whose first row is to be used as
        the new column names and then cleaned.

    Returns:
        pd.DataFrame: The modified DataFrame with cleaned column names, removed first row, and reset index.
    """

    # Create a dictionary to keep track of the count of each column name
    column_counts = {}
    # Create a list to store the new column names
    new_columns = []
    # Iterate over the columns
    for col in dataframe.iloc[0]:
        # Clean up the column name
        col = re.sub(r'\s+', ' ', str(col).strip())
        # If the column name is in the dictionary, increment the count
        if col in column_counts:
            column_counts[col] += 1
        # Otherwise, add the column name to the dictionary with a count of 1
        else:
            column_counts[col] = 1
        # If the count is greater than 1, append it to the column name
        if column_counts[col] > 1:
            col = f'{col}.{column_counts[col] - 1}'
        # Add the column name to the list of new column names
        new_columns.append(col)
    # Assign the new column names to the DataFrame
    dataframe.columns = new_columns
    # Remove the first row (which contains the original headers)
    dataframe = dataframe.iloc[1:]
    # Reset the index
    dataframe = dataframe.reset_index(drop=True)

    return dataframe


def fetch_substrate(archive, sample_id, substrate_id, logger):
    from nomad.datamodel.context import ClientContext, ServerContext
    from nomad.app.v1.models.models import User
    from nomad.search import search

    substrate_reference_str = None
    search_result = search(
        owner='all',
        query={
            'results.eln.sections:any': ['SubstrateMovpe', 'Substrate'],
            'results.eln.lab_ids:any': [substrate_id],
        },
        user_id=archive.metadata.main_author.user_id,
    )
    if not search_result.data:
        logger.warn(
            f'Substrate entry [{substrate_id}] was not found, upload and reprocess to reference it in ThinFilmStack entry [{sample_id}]'
        )
        return None
    if len(search_result.data) > 1:
        logger.warn(
            f'Found {search_result.pagination.total} entries with lab_id: '
            f'"{substrate_id}". Will use the first one found.'
        )
        return None
    if len(search_result.data) >= 1:
        upload_id = search_result.data[0]['upload_id']
        from nomad.files import UploadFiles
        from nomad.app.v1.routers.uploads import get_upload_with_read_access

        upload_files = UploadFiles.get(upload_id)

        substrate_context = ServerContext(
            get_upload_with_read_access(
                upload_id,
                User(
                    is_admin=True,
                    user_id=archive.metadata.main_author.user_id,
                ),
                include_others=True,
            )
        )

        if upload_files.raw_path_is_file(substrate_context.raw_path()):
            substrate_reference_str = f"../uploads/{search_result.data[0]['upload_id']}/archive/{search_result.data[0]['entry_id']}#data"
            return substrate_reference_str
        else:
            logger.warn(
                f"The path '../uploads/{search_result.data[0]['upload_id']}/archive/{search_result.data[0]['entry_id']}#data' is not a file, upload and reprocess to reference it in ThinFilmStack entry [{sample_id}]"
            )
            return None


def populate_sources(line_number, growth_run_file: pd.DataFrame):
    """
    Populate the Bubbler object from the growth run file
    """
    from imem_nomad_plugin.movpe.schema import BubblerSourceIMEM

    sources = []
    bubbler_quantities = [
        'Bubbler Temp',
        'Bubbler Pressure',
        'Partial Pressure',
        'Bubbler Dilution',
        'Inject',
        'Bubbler Flow',
        'Bubbler Material',
    ]
    i = 0
    while True:
        if all(
            f"{key}{'' if i == 0 else '.' + str(i)}" in growth_run_file.columns
            for key in bubbler_quantities
        ):
            sources.append(
                BubblerSourceIMEM(
                    name=growth_run_file.get(
                        f"Bubbler Material{'' if i == 0 else '.' + str(i)}", ''
                    )[line_number],
                    material=[
                        PureSubstanceComponent(
                            substance_name=growth_run_file.get(
                                f"Bubbler Material{'' if i == 0 else '.' + str(i)}", ''
                            )[line_number],
                            pure_substance=PureSubstanceSection(
                                name=growth_run_file.get(
                                    f"Bubbler Material{'' if i == 0 else '.' + str(i)}",
                                    '',
                                )[line_number]
                            ),
                        ),
                    ],
                    vapor_source=BubblerEvaporator(
                        temperature=Temperature(
                            set_value=pd.Series(
                                [
                                    growth_run_file.get(
                                        f"Bubbler Temp{'' if i == 0 else '.' + str(i)}",
                                        0,
                                    )[line_number]
                                ]
                            )
                            * ureg('celsius').to('kelvin').magnitude,
                        ),
                        pressure=Pressure(
                            set_value=pd.Series(
                                [
                                    growth_run_file.get(
                                        f"Bubbler Pressure{'' if i == 0 else '.' + str(i)}",
                                        0,
                                    )[line_number]
                                ]
                            )
                            * ureg('mbar').to('pascal').magnitude,
                        ),
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
                        total_flow_rate=VolumetricFlowRate(
                            set_value=pd.Series(
                                [
                                    growth_run_file.get(
                                        f"Bubbler Flow{'' if i == 0 else '.' + str(i)}",
                                        0,
                                    )[line_number]
                                ]
                            )
                            * ureg('cm **3 / minute')
                            .to('meter ** 3 / second')
                            .magnitude,
                        ),
                        dilution=growth_run_file.get(
                            f"Bubbler Dilution{'' if i == 0 else '.' + str(i)}", 0
                        )[line_number],
                        inject=growth_run_file.get(
                            f"Inject{'' if i == 0 else '.' + str(i)}", 0
                        )[line_number],
                    ),
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
                ),
            )

            i += 1
        else:
            break
    return sources


def populate_gas_source(line_number, growth_run_file: pd.DataFrame):
    """
    Populate the GasSource object from the growth run file
    """
    from imem_nomad_plugin.movpe.schema import GasSourceIMEM, GasLineIMEM

    gas_sources = []
    gas_source_quantities = [
        'Gas Cylinder Material',
        'Dilution in Cylinder',
        'Flow from MFC',
        'Effective  Flow',
        'Gas Partial Pressure',
        'Cylinder Pressure',
        'Gas Valve',
    ]
    i = 0
    while True:
        if all(
            f"{key}{'' if i == 0 else '.' + str(i)}" in growth_run_file.columns
            for key in gas_source_quantities
        ):
            gas_sources.append(
                GasSourceIMEM(
                    name=growth_run_file.get(
                        f"Gas Cylinder Material{'' if i == 0 else '.' + str(i)}", ''
                    )[line_number],
                    dilution_in_cylinder=growth_run_file.get(
                        f"Dilution in Cylinder{'' if i == 0 else '.' + str(i)}", ''
                    )[line_number],
                    gas_valve=growth_run_file.get(
                        f"Gas Valve{'' if i == 0 else '.' + str(i)}", ''
                    )[line_number],
                    material=[
                        PureSubstanceComponent(
                            substance_name=growth_run_file.get(
                                f"Gas Cylinder Material{'' if i == 0 else '.' + str(i)}",
                                '',
                            )[line_number],
                            pure_substance=PureSubstanceSection(
                                name=growth_run_file.get(
                                    f"Gas Cylinder Material{'' if i == 0 else '.' + str(i)}",
                                    '',
                                )[line_number]
                            ),
                        ),
                    ],
                    vapor_source=GasLineIMEM(
                        pressure=Pressure(
                            set_value=pd.Series(
                                [
                                    growth_run_file.get(
                                        f"Cylinder Pressure{'' if i == 0 else '.' + str(i)}",
                                        0,
                                    )[line_number]
                                ]
                            )
                            * ureg('mbar').to('pascal').magnitude,
                        ),
                        # precursor_partial_pressure=PartialVaporPressure(
                        #     set_value=pd.Series(
                        #         [
                        #             growth_run_file.get(
                        #                 f"Gas Partial Pressure{'' if i == 0 else '.' + str(i)}",
                        #                 0,
                        #             )[line_number]
                        #         ]
                        #     ),
                        # ),
                        total_flow_rate=VolumetricFlowRate(
                            set_value=pd.Series(
                                [
                                    growth_run_file.get(
                                        f"Flow from MFC{'' if i == 0 else '.' + str(i)}",
                                        0,
                                    )[line_number]
                                ]
                            ),
                        ),
                        effective_flow_rate=VolumetricFlowRate(
                            set_value=pd.Series(
                                [
                                    growth_run_file.get(
                                        f"Effective  Flow{'' if i == 0 else '.' + str(i)}",
                                        0,
                                    )[line_number]
                                ]
                            ),
                        ),
                    ),
                    # vapor_molar_flow_rate=MolarFlowRate(
                    #     set_value=pd.Series(
                    #         [
                    #             growth_run_file.get(
                    #                 f"Gas Molar Flux{'' if i == 0 else '.' + str(i)}",
                    #                 0,
                    #             )[line_number]
                    #         ]
                    #     )
                    #     * ureg('mol / minute').to('mol / second').magnitude,
                    # ),
                )
            )
            i += 1
        else:
            break
    return gas_sources


def populate_element(line_number, substrates_file: pd.DataFrame):
    """
    Populate the GasSource object from the growth run file
    """
    elements = []
    elements_quantities = [
        'Elements',
    ]
    i = 0
    while True:
        if all(
            f"{key}{'' if i == 0 else '.' + str(i)}" in substrates_file.columns
            for key in elements_quantities
        ):
            element = substrates_file.get(
                f"Elements{'' if i == 0 else '.' + str(i)}", ''
            )[line_number]
            if not pd.isna(element):
                elements.append(
                    ElementalComposition(
                        element=element,
                    )
                )
            i += 1
        else:
            break
    return elements


def populate_dopant(line_number, substrates_file: pd.DataFrame):
    """
    Populate the GasSource object from the growth run file
    """
    dopants = []
    dopant_quantities = [
        'Doping species',
        'Doping Level',
    ]
    i = 0
    while True:
        if all(
            f"{key}{'' if i == 0 else '.' + str(i)}" in substrates_file.columns
            for key in dopant_quantities
        ):
            doping_species = substrates_file.get(
                f"Doping species{'' if i == 0 else '.' + str(i)}", ''
            )[line_number]
            doping_level = substrates_file.get(
                f"Doping Level{'' if i == 0 else '.' + str(i)}", ''
            )[line_number]
            if not pd.isna(doping_species):
                dopants.append(
                    Dopant(
                        element=doping_species,
                        doping_level=doping_level,
                    )
                )
            i += 1
        else:
            break
    return dopants


def is_activity_section(section):
    return any('Activity' in i.label for i in section.m_def.all_base_sections)


def handle_section(section):
    if hasattr(section, 'reference') and is_activity_section(section.reference):
        return [ExperimentStep(activity=section.reference, name=section.reference.name)]
    if section.m_def.label == 'CharacterizationMovpeIMEM':
        sub_sect_list = []
        for sub_section in vars(section).values():
            if isinstance(sub_section, list):
                for item in sub_section:
                    if hasattr(item, 'reference') and is_activity_section(
                        item.reference
                    ):
                        sub_sect_list.append(
                            ExperimentStep(
                                activity=item.reference, name=item.reference.name
                            )
                        )
        return sub_sect_list
    if not hasattr(section, 'reference') and is_activity_section(section):
        return [ExperimentStep(activity=section, name=section.name)]
