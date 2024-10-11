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

from nomad.config.models.plugins import ParserEntryPoint
from pydantic import Field


class MovpeParserEntryPoint(ParserEntryPoint):
    def load(self):
        from imem_nomad_plugin.movpe.growth_excel.parser import ParserMovpeIMEM

        return ParserMovpeIMEM(**self.dict())


movpe_growth_excel_parser = MovpeParserEntryPoint(
    name='MovpeParser',
    description='Parser defined using the new plugin mechanism.',
    mainfile_name_re=r'.+\.xlsx',
    mainfile_mime_re=r'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
    mainfile_contents_dict={
        'Overview': {'__has_all_keys': ['Sample', 'Substrate T', 'VI III Ratio']},
        'Substrate': {'__has_all_keys': ['Substrates', 'Orientation']},
        'GrowthRun': {
            '__has_all_keys': ['Name', 'Flow Metal Carrier', 'Flow Oxydant Carrier']
        },
        '__comment_symbol': '#',
    },
)
