import numpy as np
import yaml
import json
import plotly.express as px
from plotly.subplots import make_subplots
from nomad.datamodel.metainfo.basesections import (
    Activity,
    ActivityStep,
    System,
    Component,
    SystemComponent,
    PureSubstance,
    Process,
    PureSubstanceComponent,
    PureSubstanceSection,
    EntityReference,
    CompositeSystemReference,
    PubChemPureSubstanceSection,
    SectionReference,
    Experiment,
    ExperimentStep,
)
from nomad.datamodel.metainfo.annotations import (
    ELNAnnotation,
    SectionProperties,
)


from nomad.parsing.tabular import TableData
from structlog.stdlib import (
    BoundLogger,
)
from nomad.datamodel.metainfo.annotations import (
    ELNAnnotation,
    ELNComponentEnum,
)

from nomad.metainfo import (
    Package,
    Quantity,
    SubSection,
    MEnum,
    Datetime,
    Section,
    Reference,
)
from nomad.datamodel.data import EntryData, ArchiveSection, Author

from nomad.datamodel.metainfo.plot import PlotSection, PlotlyFigure
from nomad.datamodel.metainfo.workflow import (
    Link,
    Task,
    Workflow,
)

from nomad_material_processing import (
    SubstrateReference,
    CrystallineSubstrate,
    Miscut,
    SubstrateCrystalProperties,
    Geometry,
    ThinFilm,
    ThinFilmStack,
    ThinFilmStackReference,
    Parallelepiped,
)
from nomad_material_processing.vapor_deposition import (
    VaporDeposition,
    VaporDepositionStep,
    SampleParameters,
    ChamberEnvironment,
    SubstrateHeater,
    Pressure,
    Temperature,
    MolarFlowRate,
    VolumetricFlowRate,
    GasFlow,
)

from nomad_material_processing.vapor_deposition.cvd import (
    BubblerEvaporator,
    FlashEvaporator,
    CVDSource,
    Rotation,
    GasSupply,
    GasLine,
)

from nomad_measurements import (
    ActivityReference,
)

from nomad_measurements.xrd import ELNXRayDiffraction

from nomad.config import config

from lakeshore_nomad_plugin.hall.schema import HallMeasurement

from imem_nomad_plugin.utils import (
    create_archive,
    handle_section,
)
from imem_nomad_plugin.general.schema import (
    IMEMMOVPECategory,
    SubstratePreparationStepReference,
    SampleCutIMEM,
)
from imem_nomad_plugin.characterization.schema import AFMmeasurement, LightMicroscope

# m_package = Package(name="movpe_IMEM")

from nomad.metainfo import (
    SchemaPackage,
)

configuration = config.get_plugin_entry_point('imem_nomad_plugin.movpe:movpe_schema')

m_package = SchemaPackage()


class Shape(Parallelepiped):
    diameter = Quantity(
        type=float,
        description='The diamater',
        a_eln=ELNAnnotation(
            component=ELNComponentEnum.NumberEditQuantity,
            defaultDisplayUnit='millimeter',
        ),
        unit='meter',
    )


class BubblerPrecursor(PureSubstance, EntryData):
    """
    A precursor already loaded in a bubbler.
    To calculate the vapor pressure the Antoine equation is used.
    log10(p) = A - [B / (T + C)]
    It is a mathematical expression (derived from the Clausius-Clapeyron equation)
    of the relation between the vapor pressure (p) and the temperature (T) of pure substances.
    """

    m_def = Section(categories=[IMEMMOVPECategory])
    name = Quantity(
        type=str,
        description='FILL',
        a_eln=ELNAnnotation(component='StringEditQuantity', label='Substance Name'),
    )
    cas_number = Quantity(
        type=str,
        description='FILL',
        a_eln=ELNAnnotation(component='StringEditQuantity', label='CAS number'),
    )
    weight = Quantity(
        type=np.float64,
        description="""
        Weight of precursor and bubbler.
        Attention: Before weighing bubblers,
        all gaskets and corresponding caps must be attached!
        """,
        a_eln=ELNAnnotation(
            component='NumberEditQuantity',
            defaultDisplayUnit='gram',
        ),
        unit='kg',
    )
    weight_difference = Quantity(
        type=np.float64,
        description='Weight when the bubbler is exhausted.',
        a_eln=ELNAnnotation(
            component='NumberEditQuantity',
            defaultDisplayUnit='gram',
        ),
        unit='kg',
    )
    total_comsumption = Quantity(
        type=np.float64,
        description='FILL DESCRIPTION.',
        a_eln=ELNAnnotation(
            component='NumberEditQuantity',
            defaultDisplayUnit='gram',
        ),
        unit='kg',
    )
    a_parameter = Quantity(
        type=np.float64,
        description='The A parameter of Antoine equation. Dimensionless.',
        a_eln=ELNAnnotation(
            component='NumberEditQuantity',
            defaultDisplayUnit='millimeter',
        ),
        unit='millimeter',
    )
    b_parameter = Quantity(
        type=np.float64,
        description='The B parameter of Antoine equation. Temperature units.',
        a_eln=ELNAnnotation(
            component='NumberEditQuantity',
            defaultDisplayUnit='celsius',
        ),
        unit='kelvin',
    )
    c_parameter = Quantity(
        type=np.float64,
        description='The C parameter of Antoine equation. Temperature units.',
        a_eln=ELNAnnotation(
            component='NumberEditQuantity',
            defaultDisplayUnit='celsius',
        ),
        unit='kelvin',
    )
    information_sheet = Quantity(
        type=str,
        description='pdf files containing certificate and other documentation',
        a_browser={'adaptor': 'RawFileAdaptor'},
        a_eln=ELNAnnotation(
            component='FileEditQuantity',
        ),
    )


class Cylinder(Geometry):
    """
    Class autogenerated from yaml schema.
    """

    m_def = Section()
    height = Quantity(
        type=np.float64,
        description='docs',
        a_eln=ELNAnnotation(
            component='NumberEditQuantity',
            defaultDisplayUnit='nanometer',
        ),
        unit='nanometer',
    )
    radius = Quantity(
        type=np.float64,
        description='docs',
        a_eln=ELNAnnotation(
            component='NumberEditQuantity',
            defaultDisplayUnit='millimeter',
        ),
        unit='millimeter',
    )
    lower_cap_radius = Quantity(
        type=np.float64,
        description='docs',
        a_eln=ELNAnnotation(
            component='NumberEditQuantity',
            defaultDisplayUnit='millimeter',
        ),
        unit='millimeter',
    )
    upper_cap_radius = Quantity(
        type=np.float64,
        description='docs',
        a_eln=ELNAnnotation(
            component='NumberEditQuantity',
            defaultDisplayUnit='millimeter',
        ),
        unit='millimeter',
    )
    cap_surface_area = Quantity(
        type=np.float64,
        description='docs',
        a_eln=ELNAnnotation(
            component='NumberEditQuantity',
            defaultDisplayUnit='millimeter ** 2',
        ),
        unit='millimeter ** 2',
    )
    lateral_surface_area = Quantity(
        type=np.float64,
        description='docs',
        a_eln=ELNAnnotation(
            component='NumberEditQuantity',
            defaultDisplayUnit='millimeter ** 2',
        ),
        unit='millimeter ** 2',
    )


class MiscutMovpe(Miscut):
    """
    The miscut in a crystalline substrate refers to
    the intentional deviation from a specific crystallographic orientation,
    commonly expressed as the angular displacement of a crystal plane.
    """

    m_def = Section(label='Miscut')

    b_angle = Quantity(
        type=float,
        description='crystallographic orientation of the substrate in [hkl]',
        a_eln=ELNAnnotation(
            component='NumberEditQuantity',
        ),
        a_tabular={
            'name': 'Substrate/Miscut b angle',
            # "unit": "deg"
        },
        unit='deg',
    )
    angle = Quantity(
        type=float,
        description='angular displacement from crystallographic orientation of the substrate',
        a_eln=ELNAnnotation(
            component='NumberEditQuantity',
            defaultDisplayUnit='deg',
            label='c angle',
        ),
        unit='deg',
        a_tabular={
            'name': 'Substrate/Miscut c angle',
            # "unit": "deg"
        },
    )
    angle_deviation = Quantity(
        type=float,
        description='uncertainty on the angular displacement',
        a_eln=ELNAnnotation(
            component='NumberEditQuantity',
            defaultDisplayUnit='deg',
            label='c angle deviation',
        ),
        unit='deg',
    )
    orientation = Quantity(
        type=str,
        description='crystallographic orientation of the substrate in [hkl]',
        a_eln=ELNAnnotation(
            component='StringEditQuantity',
        ),
        a_tabular={'name': 'Substrate/Miscut c Orientation'},
    )


class SubstrateCrystalPropertiesMovpe(SubstrateCrystalProperties):
    """
    Characteristics arising from the ordered arrangement of atoms in a crystalline structure.
    These properties are defined by factors such as crystal symmetry, lattice parameters,
    and the specific arrangement of atoms within the crystal lattice.
    """

    m_def = Section(label='CrystalProperties')
    orientation = Quantity(
        type=str,
        a_eln=ELNAnnotation(
            component='StringEditQuantity',
        ),
        a_tabular={'name': 'Substrate/Orientation'},
    )
    miscut = SubSection(section_def=MiscutMovpe)


class SubstrateMovpe(CrystallineSubstrate, EntryData):
    """
    Class autogenerated from yaml schema.
    """

    m_def = Section(
        label_quantity='lab_id', categories=[IMEMMOVPECategory], label='Substrate'
    )
    tags = Quantity(
        type=str,
        description='FILL',
        a_eln=ELNAnnotation(
            component='StringEditQuantity',
            label='Box ID',
        ),
        a_tabular={'name': 'Substrate/Substrate Box'},
    )
    annealing = Quantity(
        type=bool,
        description='Annealing',
        a_eln=ELNAnnotation(
            component='BoolEditQuantity',
        ),
    )
    cleaning = Quantity(
        type=bool,
        description='Cleaning',
        a_eln=ELNAnnotation(
            component='BoolEditQuantity',
        ),
    )
    regrowth = Quantity(
        type=bool,
        description='Regrowth',
        a_eln=ELNAnnotation(
            component='BoolEditQuantity',
        ),
    )
    information_sheet = Quantity(
        type=str,
        description='pdf files containing certificate and other documentation',
        a_browser={'adaptor': 'RawFileAdaptor'},
        a_eln=ELNAnnotation(
            component='FileEditQuantity',
        ),
    )
    description = Quantity(
        type=str,
        description='description',
        a_eln=ELNAnnotation(
            component='StringEditQuantity',
            label='Notes',
        ),
    )


class ThinFilmMovpe(ThinFilm, EntryData):
    """
    Class autogenerated from yaml schema.
    """

    m_def = Section(
        label_quantity='lab_id',
        categories=[IMEMMOVPECategory],
        label='ThinFilmMovpe',
    )
    lab_id = Quantity(
        type=str,
        description='the Sample created in the current growth',
        a_tabular={'name': 'GrowthRun/Sample Name'},
        a_eln=ELNAnnotation(
            component='StringEditQuantity',
            label='Grown Sample ID',
        ),
    )
    test_quantities = Quantity(
        type=str,
        description='Test quantity',
        a_eln=ELNAnnotation(
            component='StringEditQuantity',
        ),
    )


class ThinFilmStackMovpeIMEM(ThinFilmStack, EntryData):
    """
    Class autogenerated from yaml schema.
    """

    m_def = Section(
        label_quantity='lab_id',
        categories=[IMEMMOVPECategory],
        label='ThinFilmStackMovpe',
    )
    lab_id = Quantity(
        type=str,
        description='the Sample created in the current growth',
        a_tabular={'name': 'GrowthRun/Sample Name'},
        a_eln=ELNAnnotation(
            component='StringEditQuantity',
            label='Grown Sample ID',
        ),
    )
    parent_sample = SubSection(
        description="""
        the parent sample of the current sample.
        """,
        section_def=ThinFilmStackReference,
    )


class ThinFilmStackMovpeReference(ThinFilmStackReference):
    """
    A section used for referencing a Grown Sample.
    """

    lab_id = Quantity(
        type=str,
        description='the Sample created in the current growth',
        a_tabular={'name': 'GrowthRun/Sample Name'},
        a_eln=ELNAnnotation(
            component='StringEditQuantity',
            label='Grown Sample ID',
        ),
    )
    reference = Quantity(
        type=ThinFilmStackMovpeIMEM,
        description='A reference to a NOMAD `ThinFilmStackMovpe` entry.',
        a_eln=ELNAnnotation(
            component='ReferenceEditQuantity',
            label='ThinFilmStackMovpe Reference',
        ),
    )

    def normalize(self, archive, logger: BoundLogger) -> None:
        """
        The normalizer for the `ThinFilmStackMovpeReference` class.
        """
        super(ThinFilmStackMovpeReference, self).normalize(archive, logger)


class SystemComponentIMEM(SystemComponent):
    """
    A section for describing a system component and its role in a composite system.
    """

    molar_concentration = Quantity(
        type=np.float64,
        description='The solvent for the current substance.',
        unit='mol/liter',
        a_eln=dict(component='NumberEditQuantity', defaultDisplayUnit='mol/liter'),
        a_tabular={
            'name': 'Precursors/Molar conc',
            # "unit": "gram"
        },
    )
    system = Quantity(
        type=Reference(System.m_def),
        description='A reference to the component system.',
        a_eln=dict(component='ReferenceEditQuantity'),
    )


class PrecursorsPreparationIMEM(Process, EntryData):
    """
    Class autogenerated from yaml schema.
    """

    m_def = Section(
        a_eln={
            'hide': [
                'instruments',
                'steps',
                'samples',
            ]
        },
        label_quantity='name',
        categories=[IMEMMOVPECategory],
        label='PrecursorsPreparation',
    )
    data_file = Quantity(
        type=str,
        description='Upload here the spreadsheet file containing the deposition control data',
        a_browser={'adaptor': 'RawFileAdaptor'},
        a_eln={'component': 'FileEditQuantity'},
    )
    lab_id = Quantity(
        type=str,
        description='FILL',
        a_tabular={'name': 'Precursors/Sample ID'},
        a_eln={'component': 'StringEditQuantity', 'label': 'Sample ID'},
    )
    name = Quantity(
        type=str,
        description='FILL',
        a_tabular={'name': 'Precursors/number'},
        a_eln={
            'component': 'StringEditQuantity',
        },
    )
    description = Quantity(
        type=str,
        a_eln={'component': 'StringEditQuantity'},
    )
    flow_titanium = Quantity(  # TODO make this a single flow
        type=np.float64,
        description='FILL THE DESCRIPTION',
        a_tabular={'name': 'Precursors/Set flow Ti'},
        a_eln={'component': 'NumberEditQuantity', 'defaultDisplayUnit': 'ml / minute'},
        unit='ml / minute',
    )
    flow_calcium = Quantity(
        type=np.float64,
        description='FILL THE DESCRIPTION',
        a_tabular={'name': 'Precursors/Set flow Ca'},
        a_eln={'component': 'NumberEditQuantity', 'defaultDisplayUnit': 'ml / minute'},
        unit='ml / minute',
    )
    # precursors = SubSection(
    #     section_def=SystemComponent,
    #     description="""
    #     A precursor used in MOVPE. It can be a solution, a gas, or a solid.
    #     """,
    #     repeats=True,
    # )
    components = SubSection(
        description="""
        A list of all the components of the composite system containing a name, reference
        to the system section and mass of that component.
        """,
        section_def=Component,
        repeats=True,
    )


class PrecursorsPreparationIMEMReference(ActivityReference):
    """
    A section used for referencing a PrecursorsPreparationIMEM.
    """

    m_def = Section(
        label='PrecursorsPreparationReference',
    )
    reference = Quantity(
        type=PrecursorsPreparationIMEM,
        description='A reference to a NOMAD `PrecursorsPreparationIMEM` entry.',
        a_eln=ELNAnnotation(
            component='ReferenceEditQuantity',
            label='PrecursorsPreparationIMEM Reference',
        ),
    )


class InSituMonitoringReference(SectionReference):
    """
    A section used for referencing a InSituMonitoring.
    """

    reference = Quantity(
        type=ArchiveSection,
        description='A reference to a NOMAD `InSituMonitoring` entry.',
        a_eln=ELNAnnotation(
            component='ReferenceEditQuantity',
            label='In situ Monitoring Reference',
        ),
    )


class HallMeasurementReference(SectionReference):
    """
    A section used for referencing a HallMeasurement.
    The class is taken from the dedicated Lakeshore plugin
    """

    reference = Quantity(
        type=ArchiveSection,
        description='A reference to a NOMAD `HallMeasurement` entry.',
        a_eln=ELNAnnotation(
            component='ReferenceEditQuantity',
            label='Hall Measurement Reference',
        ),
    )


class SubstrateMovpeReference(SubstrateReference):
    """
    A section for describing a system component and its role in a composite system.
    """

    lab_id = Quantity(
        type=str,
        a_eln=ELNAnnotation(
            component=ELNComponentEnum.StringEditQuantity,
            label='Substrate ID',
        ),
    )
    reference = Quantity(
        type=SubstrateMovpe,
        a_eln=ELNAnnotation(
            component=ELNComponentEnum.ReferenceEditQuantity,
            label='Substrate',
        ),
    )


class SubstrateInventory(EntryData, TableData):
    """
    Class autogenerated from yaml schema.
    """

    m_def = Section(
        a_eln=None,
        categories=[IMEMMOVPECategory],
        label='SubstrateInventory',
    )
    data_file = Quantity(
        type=str,
        description='Upload here the spreadsheet file containing the substrates data',
        # a_tabular_parser={
        #     "parsing_options": {"comment": "#"},
        #     "mapping_options": [
        #         {
        #             "mapping_mode": "row",
        #             "file_mode": "multiple_new_entries",
        #             "sections": ["substrates"],
        #         }
        #     ],
        # },
        a_browser={'adaptor': 'RawFileAdaptor'},
        a_eln={'component': 'FileEditQuantity'},
    )
    substrates = SubSection(
        section_def=SubstrateMovpeReference,
        repeats=True,
    )
    steps = SubSection(
        section_def=SubstratePreparationStepReference,
        repeats=True,
    )


class AFMmeasurementReference(SectionReference):
    """
    A section used for referencing a AFMmeasurement.
    """

    reference = Quantity(
        type=AFMmeasurement,
        description='A reference to a NOMAD `AFMmeasurement` entry.',
        a_eln=ELNAnnotation(
            component='ReferenceEditQuantity',
            label='AFM Measurement Reference',
        ),
    )


class LiMimeasurementReference(SectionReference):
    """
    A section used for referencing a LightMicroscope.
    """

    reference = Quantity(
        type=LightMicroscope,
        description='A reference to a NOMAD `LightMicroscope` entry.',
        a_eln=ELNAnnotation(
            component='ReferenceEditQuantity',
            label='Light Microscope Measurement Reference',
        ),
    )


class XRDmeasurementReference(SectionReference):
    """
    A section used for referencing a LightMicroscope.
    """

    sample_id = Quantity(
        type=str,
        description='The sample to be linked within the XRD measurement',
        a_eln=ELNAnnotation(
            component=ELNComponentEnum.StringEditQuantity,
        ),
    )
    reference = Quantity(
        type=ELNXRayDiffraction,
        description='A reference to a NOMAD `ELNXRayDiffraction` entry.',
        a_eln=ELNAnnotation(
            component='ReferenceEditQuantity',
            label='XRD Measurement Reference',
        ),
    )
    phase = Quantity(
        type=str,
        description='Phase type obtained from HRXRD',
        a_eln=ELNAnnotation(
            component=ELNComponentEnum.StringEditQuantity,
        ),
    )
    peak_position_2theta = Quantity(
        type=np.float64,
        description='Peak Position - 2theta',
        a_eln=ELNAnnotation(
            component=ELNComponentEnum.NumberEditQuantity,
        ),
        unit='degree',
    )
    peak_fwhm_2theta = Quantity(
        type=np.float64,
        description='Peak Position - 2theta',
        a_eln=ELNAnnotation(
            component=ELNComponentEnum.NumberEditQuantity,
        ),
        unit='degree',
    )
    peak_position_omega = Quantity(
        type=np.float64,
        description='Peak Position - Omega',
        a_eln=ELNAnnotation(
            component=ELNComponentEnum.NumberEditQuantity,
        ),
        unit='degree',
    )
    peak_fwhm_rocking_curve = Quantity(
        type=str,
        description='Peak FWHM Rocking Curve',
        a_eln=ELNAnnotation(
            component=ELNComponentEnum.StringEditQuantity,
        ),
    )
    reflection = Quantity(
        type=str,
        description='Peak FWHM Rocking Curve',
        a_eln={'component': 'StringEditQuantity'},
    )
    description = Quantity(
        type=str,
        description='Notes and comments.',
        a_eln=ELNAnnotation(
            component=ELNComponentEnum.RichTextEditQuantity,
        ),
    )

    def normalize(self, archive, logger):
        super(XRDmeasurementReference, self).normalize(archive, logger)
        if (
            hasattr(self, 'reference')
            and self.reference is not None
            and hasattr(self, 'sample_id')
        ):
            from nomad.datamodel.context import ServerContext
            from nomad.app.v1.routers.uploads import get_upload_with_read_access
            from nomad.datamodel.data import User

            # xrd_context = ServerContext(
            #     get_upload_with_read_access(
            #         archive.m_context.upload_id,
            #         User(
            #             is_admin=True,
            #             user_id=archive.metadata.main_author.user_id,
            #         ),
            #         include_others=True,
            #     )
            # )

            with archive.m_context.raw_file(
                self.reference.m_parent.metadata.mainfile, 'r'
            ) as xrd_file:
                updated_xrd_file = json.load(xrd_file)
                updated_xrd_file['data']['samples'] = [
                    CompositeSystemReference(
                        lab_id=self.sample_id,
                    ).m_to_dict()
                ]

            create_archive(
                updated_xrd_file,
                archive.m_context,
                self.reference.m_parent.metadata.mainfile,
                'json',
                logger,
                overwrite=True,
            )


class CharacterizationMovpeIMEM(ArchiveSection):
    """
    A wrapped class to gather all the characterization methods in MOVPE
    """

    xrd = SubSection(
        section_def=XRDmeasurementReference,
        repeats=True,
    )
    hall = SubSection(
        section_def=HallMeasurementReference,
        repeats=True,
    )
    afm = SubSection(
        section_def=AFMmeasurementReference,
        repeats=True,
    )
    light_microscopy = SubSection(
        section_def=LiMimeasurementReference,
        repeats=True,
    )


class ShaftTemperature(Temperature):
    """
    Central shaft temperature (to hold the susceptor)
    """

    pass


class FilamentTemperature(Temperature):
    """
    heating filament temperature
    """

    pass


class LayTecTemperature(Temperature):
    """
    Central shaft temperature (to hold the susceptor)
    """

    pass


class BubblerSourceIMEM(CVDSource):
    vapor_source = SubSection(
        section_def=BubblerEvaporator,
    )


class FlashSourceIMEM(CVDSource):
    vapor_source = SubSection(
        section_def=FlashEvaporator,
        description="""
        Example: A heater, a filament, a laser, a bubbler, etc.
        """,
    )


class GasLineIMEM(GasLine):
    """
    A gas line used in MOVPE
    """

    effective_flow_rate = SubSection(
        section_def=VolumetricFlowRate,
        description="""
        Effective flow rate, to be defined better.
        """,
    )


class GasSourceIMEM(CVDSource):
    dilution_in_cylinder = Quantity(
        type=np.float64,
        description='FILL THE DESCRIPTION',
        a_eln={'component': 'NumberEditQuantity'},
    )
    gas_valve = Quantity(
        type=bool,
        description='is the valve open?',
        a_eln=ELNAnnotation(
            component='BoolEditQuantity',
        ),
    )
    vapor_source = SubSection(
        section_def=GasLineIMEM,
    )


class GasFlowMovpe(GasFlow):
    gas = SubSection(
        section_def=PureSubstanceSection,
    )
    flow_rate = SubSection(
        section_def=VolumetricFlowRate,
        label='Push Flow Rate',
    )
    purge_flow_rate = SubSection(
        section_def=VolumetricFlowRate,
    )


class ChamberEnvironmentMovpe(ChamberEnvironment):
    uniform_gas_flow_rate = SubSection(
        section_def=VolumetricFlowRate,
    )
    pressure = SubSection(
        section_def=Pressure,
    )
    throttle_valve = SubSection(
        section_def=Pressure,
    )
    rotation = SubSection(
        section_def=Rotation,
    )
    heater = SubSection(
        section_def=SubstrateHeater,
    )


class SampleParametersMovpe(SampleParameters):
    m_def = Section(
        a_eln=ELNAnnotation(
            properties=SectionProperties(
                order=[
                    'shaft_temperature',
                    'filament_temperature',
                    'laytec_temperature',
                    'substrate_temperature',
                    'in_situ_reflectance',
                    'growth_rate',
                    'layer',
                    'substrate',
                ],
            ),
        ),
        a_plotly_graph_object=[
            {
                'label': 'filament temperature',
                'index': 1,
                'dragmode': 'pan',
                'data': {
                    'type': 'scattergl',
                    'line': {'width': 2},
                    'marker': {'size': 6},
                    'mode': 'lines+markers',
                    'name': 'Filament Temperature',
                    'x': '#filament_temperature/time',
                    'y': '#filament_temperature/value',
                },
                'layout': {
                    'title': {'text': 'Filament Temperature'},
                    'xaxis': {
                        'showticklabels': True,
                        'fixedrange': True,
                        'ticks': '',
                        'title': {'text': 'Process time [min]'},
                        # "showline": True,
                        'linewidth': 1,
                        'linecolor': 'black',
                        'mirror': True,
                    },
                    'yaxis': {
                        'showticklabels': True,
                        'fixedrange': True,
                        'ticks': '',
                        'title': {'text': 'Temperature [°C]'},
                        # "showline": True,
                        'linewidth': 1,
                        'linecolor': 'black',
                        'mirror': True,
                    },
                    'showlegend': False,
                },
                'config': {
                    'displayModeBar': False,
                    'scrollZoom': False,
                    'responsive': False,
                    'displaylogo': False,
                    'dragmode': False,
                },
            },
        ],
    )
    name = Quantity(
        type=str,
        description="""
        Sample name.
        """,
        a_eln=ELNAnnotation(
            component='StringEditQuantity',
        ),
    )
    distance_to_source = Quantity(
        type=float,
        unit='meter',
        a_eln={'component': 'NumberEditQuantity', 'defaultDisplayUnit': 'millimeter'},
        description="""
        The distance between the substrate and the source.
        It is an array because multiple sources can be used.
        """,
        shape=[1],
    )
    filament_temperature = SubSection(
        section_def=FilamentTemperature,
    )
    in_situ_reflectance = SubSection(
        section_def=InSituMonitoringReference,
    )


class GrowthStepMovpeIMEM(VaporDepositionStep, PlotSection):
    """
    Growth step for MOVPE IMEM
    """

    m_def = Section(
        label='Growth Step Movpe 2',
        a_eln=None,
    )
    # name
    # step_index
    # creates_new_thin_film
    # duration
    # sources
    # sample_parameters
    # environment
    # description

    name = Quantity(
        type=str,
        description="""
        A short and descriptive name for this step.
        """,
        a_eln=ELNAnnotation(
            component='StringEditQuantity',
            label='Step name',
        ),
    )
    step_index = Quantity(
        type=str,
        description='the ID from RTG',
        a_eln={
            'component': 'StringEditQuantity',
        },
    )
    # duration = VaporDepositionStep.duration.m_copy()

    duration = Quantity(
        type=float,
        unit='second',
        a_eln=ELNAnnotation(
            component='NumberEditQuantity',
        ),
    )

    comment = Quantity(
        type=str,
        description='description',
        a_eln={'component': 'StringEditQuantity'},
        label='Notes',
    )
    sample_parameters = SubSection(
        section_def=SampleParametersMovpe,
        repeats=True,
    )
    sources = SubSection(
        section_def=CVDSource,
        repeats=True,
    )
    environment = SubSection(
        section_def=ChamberEnvironmentMovpe,
    )
    in_situ_reflectance = SubSection(
        section_def=InSituMonitoringReference,
    )


class GrowthMovpeIMEM(VaporDeposition, EntryData):
    """
    Class autogenerated from yaml schema.
    """

    m_def = Section(
        a_eln=ELNAnnotation(
            properties=SectionProperties(
                order=[
                    'name',
                    'method',
                    'data_file',
                    'datetime',
                    'end_time',
                    'duration',
                ],
            ),
            # hide=[
            #     "instruments",
            #     "steps",
            #     "samples",
            #     "description",
            #     "location",
            #     "lab_id",
            # ],
        ),
        label_quantity='lab_id',
        categories=[IMEMMOVPECategory],
        label='Growth Process',
    )

    # datetime
    # name
    # description
    # lab_id
    # method
    method = Quantity(
        type=str,
        default='MOVPE IMEM',
    )
    data_file = Quantity(
        type=str,
        description='Upload here the spreadsheet file containing the deposition control data',
        # a_tabular_parser={
        #     "parsing_options": {"comment": "#"},
        #     "mapping_options": [
        #         {
        #             "mapping_mode": "row",
        #             "file_mode": "multiple_new_entries",
        #             "sections": ["#root"],
        #         }
        #     ],
        # },
        a_browser={'adaptor': 'RawFileAdaptor'},
        a_eln={'component': 'FileEditQuantity'},
    )
    description = Quantity(
        type=str,
        description='description',
        a_eln={'component': 'StringEditQuantity'},
        label='Notes',
    )
    recipe_id = Quantity(
        type=str,
        description='the ID from RTG',
        a_tabular={'name': 'GrowthRun/Recipe Name'},
        a_eln={'component': 'StringEditQuantity', 'label': 'Recipe ID'},
    )
    susceptor = Quantity(
        type=str,
        description="""
        material of the susceptor adaptor
        """,
        a_eln=ELNAnnotation(
            component='StringEditQuantity',
        ),
    )
    mask = Quantity(
        type=str,
        description="""
        type and size of growth map
        """,
        a_eln=ELNAnnotation(
            component='StringEditQuantity',
        ),
    )
    pocket = Quantity(
        type=str,
        description="""
        position in the growth mask
        """,
        a_eln=ELNAnnotation(
            component='StringEditQuantity',
        ),
    )
    steps = SubSection(
        section_def=GrowthStepMovpeIMEM,
        repeats=True,
    )

    def normalize(self, archive, logger):
        # for sample in self.samples:
        #     sample.normalize(archive, logger)
        # for parent_sample in self.parent_sample:
        #     parent_sample.normalize(archive, logger)
        # for substrate in self.substrate:
        #     substrate.normalize(archive, logger)

        archive.workflow2 = None
        super(GrowthMovpeIMEM, self).normalize(archive, logger)
        if self.steps is not None:
            inputs = []
            outputs = []
            for step in self.steps:
                if step.sample_parameters is not None:
                    for sample in step.sample_parameters:
                        if sample.layer is not None:
                            outputs.append(
                                Link(
                                    name=f'{sample.layer.name}',
                                    section=sample.layer.reference,
                                )
                            )
                        if sample.substrate is not None:
                            outputs.append(
                                Link(
                                    name=f'{sample.substrate.name}',
                                    section=sample.substrate.reference,
                                )
                            )
                        if (
                            sample.substrate is not None
                            and sample.substrate.reference is not None
                        ):
                            if hasattr(
                                getattr(sample.substrate.reference, 'substrate'),
                                'name',
                            ):
                                # sample.substrate.reference.substrate.reference is not None:
                                inputs.append(
                                    Link(
                                        name=f'{sample.substrate.reference.substrate.name}',
                                        section=getattr(
                                            sample.substrate.reference.substrate,
                                            'reference',
                                            None,
                                        ),
                                    )
                                )
            archive.workflow2.outputs.extend(set(outputs))
            archive.workflow2.inputs.extend(set(inputs))


class GrowthMovpeIMEMReference(ActivityReference):
    """
    A section used for referencing a GrowthMovpeIMEM.
    """

    m_def = Section(
        label='SampleCutReference',
    )
    reference = Quantity(
        type=SampleCutIMEM,
        description='A reference to a NOMAD `SampleCutIMEM` entry.',
        a_eln=ELNAnnotation(
            component='ReferenceEditQuantity',
        ),
    )


class SampleCutIMEMReference(ActivityReference):
    """
    A section used for referencing a GrowthMovpeIMEM.
    """

    m_def = Section(
        label='GrowthProcessReference',
    )
    reference = Quantity(
        type=GrowthMovpeIMEM,
        description='A reference to a NOMAD `GrowthMovpeIMEM` entry.',
        a_eln=ELNAnnotation(
            component='ReferenceEditQuantity',
        ),
    )


class ExperimentMovpeIMEM(Experiment, EntryData):
    """
    Class autogenerated from yaml schema.
    """

    m_def = Section(
        # a_eln={"hide": ["steps"]},
        categories=[IMEMMOVPECategory],
        label='MOVPE Experiment',
    )
    # lab_id
    method = Quantity(
        type=str,
    )
    data_file = Quantity(
        type=str,
        description='Upload here the spreadsheet file containing the growth data',
        a_browser={'adaptor': 'RawFileAdaptor'},
        a_eln={'component': 'FileEditQuantity'},
    )
    description = Quantity(
        type=str,
        description='description',
        a_eln=ELNAnnotation(
            component='StringEditQuantity',
            label='Notes',
        ),
    )
    substrate_temperature = Quantity(
        type=np.float64,
        description='FILL THE DESCRIPTION',
        a_eln=ELNAnnotation(
            component='NumberEditQuantity',
            defaultDisplayUnit='celsius',
        ),
        unit='kelvin',
    )
    oxygen_argon_ratio = Quantity(
        type=str,
        description='FILL THE DESCRIPTION',
        a_eln=ELNAnnotation(
            component='StringEditQuantity',
        ),
    )
    composition = Quantity(
        type=str,
        description='FILL THE DESCRIPTION',
        a_eln={
            'component': 'StringEditQuantity',
        },
    )
    precursors_preparation = SubSection(
        section_def=PrecursorsPreparationIMEMReference,
    )

    pregrowth = SubSection(
        section_def=GrowthMovpeIMEMReference,
    )
    growth_run = SubSection(
        section_def=GrowthMovpeIMEMReference,
    )
    sample_cut = SubSection(
        section_def=GrowthMovpeIMEMReference,
    )
    characterization = SubSection(section_def=CharacterizationMovpeIMEM)

    steps = SubSection(
        section_def=ActivityReference,
        repeats=True,
    )
    # growth_run_constant_parameters = SubSection(
    #     section_def=GrowthMovpe1IMEMConstantParametersReference
    # )

    def normalize(self, archive, logger):
        archive_sections = (
            attr for attr in vars(self).values() if isinstance(attr, ArchiveSection)
        )
        step_list = []
        for section in archive_sections:
            try:
                step_list.extend(handle_section(section))
            except (AttributeError, TypeError, NameError) as e:
                print(f'An error occurred: {e}')
        self.steps = [step for step in step_list if step is not None]

        archive.workflow2 = None
        super(ExperimentMovpeIMEM, self).normalize(archive, logger)

        # search_result = search(
        #     owner="user",
        #     query={
        #         "results.eln.sections:any": ["GrowthMovpe1IMEMConstantParameters"],
        #         "upload_id:any": [archive.m_context.upload_id],
        #     },
        #     pagination=MetadataPagination(page_size=10000),
        #     user_id=archive.metadata.main_author.user_id,
        # )
        # # checking if all entries are properly indexed
        # if getattr(
        #     getattr(self, "growth_run_constant_parameters", None), "lab_id", None
        # ) and not getattr(
        #     getattr(self, "growth_run_constant_parameters", None), "reference", None
        # ):
        #     found_id = False
        #     for growth_entry in search_result.data:
        #         if (
        #             self.growth_run_constant_parameters.lab_id
        #             == growth_entry["results"]["eln"]["lab_ids"][0]
        #         ):
        #             found_id = True
        #             self.growth_run_constant_parameters = GrowthMovpe1IMEMConstantParametersReference(
        #                 reference=f"../uploads/{archive.m_context.upload_id}/archive/{growth_entry['entry_id']}#data"
        #             )
        #         for search_quantities in growth_entry["search_quantities"]:
        #             if (
        #                 search_quantities["path_archive"]
        #                 == "data.substrate_temperature"
        #             ):
        #                 self.substrate_temperature = search_quantities["float_value"]
        #             if search_quantities["path_archive"] == "data.oxygen_argon_ratio":
        #                 self.oxygen_argon_ratio = search_quantities["float_value"]
        #             if search_quantities["path_archive"] == "data.composition":
        #                 self.composition = search_quantities["str_value"][0]
        #     if not found_id:
        #         logger.warning(
        #             f"The lab_id '{self.growth_run_constant_parameters.lab_id}' was not found in any 'GrowthMovpe1IMEMConstantParameters' entry in Nomad. Check if it exist and try to reference it manually."
        #         )
        # else:
        #     logger.warning(
        #         "No lab_id for 'GrowthMovpe1IMEMConstantParameters' found. The archive couldn't be referenced."
        #     )

    # def normalize(self, archive, logger: BoundLogger) -> None:
    #     '''
    #     The normalizer for the `MovpeBinaryOxidesIMEMExperiment` class.
    #     '''
    #     super(MovpeBinaryOxidesIMEMExperiment, self).normalize(archive, logger)
    ## Potential weak code in next lines:
    ## I want to get back to GrowthRun entry (already created by tabular parser)
    ## and set the "reference" quantity in grwon_samples.
    ## Here two example codes by Theodore Chang, first touches the raw file, second touches the processed file.
    #### ONE
    ## 1. get the file name of archive/entry containing grown_sample_ref
    ## 2. overwrite yaml for this entry
    ## 3. reprocess
    # grown_sample_ref.reference = f'../uploads/{archive.m_context.upload_id}/archive/{hash(archive.m_context.upload_id, filename)}#data'
    # grown_sample_archive = grown_sample_ref
    # while not isinstance(grown_sample_archive, EntryArchive):
    #     grown_sample_archive=grown_sample_archive.m_parent
    # grown_sample_file_name:str = grown_sample_archive.metadata.mainfile
    # create_archive(
    #     grown_sample_archive.m_to_dict(), archive.m_context, grown_sample_file_name, filetype, logger,bypass_check=True)
    #### TWO
    ## alternatively directly overwite the processed msg file
    # grown_sample_upload_id:str = grown_sample_archive.metadata.upload_id
    # grown_sample_entry_id:str = grown_sample_archive.metadata.entry_id
    # StagingUploadFiles(grown_sample_upload_id).write_archive(grown_sample_entry_id, grown_sample_archive.m_to_dict())


m_package.__init_metainfo__()
