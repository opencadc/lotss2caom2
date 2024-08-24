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
#  Revision: 4
#
# ***********************************************************************
#

import numpy as np
from astropy.table import Table
from caom2 import Algorithm, SimpleObservation


def _search_id_list_mock(ignore_query):
    # the query return value is pyvo.dal.tap.TAPResults, but for the purposes of testing
    # astropy.table.Table has the same behaviour
    return Table.read('mosaic_id\nP000+23\nP000+31\nP000+36\n'.split('\n'), format='csv')


# 'https://vo.astron.nl/getproduct/LoTSS-DR2/P000%2B23\t\t\timage/fits\t639106560\t0.03125\t23.3953\tP000+23_mosaic-blanked.fits\tLOFAR.HBA\t58452.6351273148\t2\tmasked_array(data=[8938, 8938]\t masked_array(data=[0.0004166669968981296,0.0004166669968981296]\tICRS\t2000.0\tSIN\tmasked_array(data=[4469.0, 4469.0]\tmasked_array(data=[0.03125, 23.3953]\tmasked_array(data=[-0.0004166669968981296, 0.0, 0.0, 0.0004166669968981296]\t120-168MHz\tm\t2.08189\t2.49827\t1.78448\tZ\tmasked_array(data=[2.032831406303902, 21.520385451900815, 2.08992790594909, 25.244422592409883, 357.9721112359749, 25.24441664211163, 358.02922050953805, 21.520379666564615]\tP000+23\thttps://vo.astron.nl/lotss_dr2/q/dlmosaic/dlmeta?ID=ivo%3A//astron.nl/%7E%3FLoTSS-DR2/P000%2B23\tmasked_array(data=[689778]\t21.12136/ede3eff0-266a-465a-af68-6f292391485f\t8eefe6e9\n'.split('\n'),
def _search_mosaic_id_mock(query_string):
    # the query return value is pyvo.dal.tap.TAPResults, and the individual masked arrays are numpy.ma.core.MaskedArray
    # but for the purposes of testing a list of dicts has the same behaviour
    # service = vo.dal.TAPService('https://vo.astron.nl/tap')
    # results = service.search(f"SELECT * from lotss_dr2.mosaics WHERE mosaic_id = '{mosaic_id}'")
    if 'P000+23' in query_string:
        result = [
            {
                'accref': 'https://vo.astron.nl/getproduct/LoTSS-DR2/P000%2B23',
                'owner': '',
                'embargo': '',
                'mime': 'image/fits',
                'accsize': np.int64(639106560),
                'centeralpha': 0.03125,
                'centerdelta': 23.3953,
                'imageTitle': 'P000+23_mosaic-blanked.fits',
                'instid': 'LOFAR.HBA',
                'dateobs': 58452.6351273148,
                'nAxes': 2,
                'pixelsize': np.array([8938, 8938]),
                'pixelScale': np.array([0.0004166669968981296, 0.0004166669968981296]),
                'refframe': 'ICRS',
                'wcs_equinox': 2000.0,
                'wcs_projection': 'SIN',
                'wcs_refpixel': np.array([4469.0, 4469.0]),
                'wcs_refvalues': np.array([0.03125, 23.3953]),
                'wcs_cdmatrix': np.array([-0.0004166669968981296, 0.0, 0.0, 0.0004166669968981296]),
                'bandpassid': '120-168MHz',
                'bandpassunit': 'm',
                'bandpassrefval': 2.08189,
                'bandpasshi': 2.49827,
                'bandpasslo': 1.78448,
                'pixflags': 'Z',
                'coverage': np.array(
                    [
                        2.032831406303902,
                        21.520385451900815,
                        2.08992790594909,
                        25.244422592409883,
                        357.9721112359749,
                        25.24441664211163,
                        358.02922050953805,
                        21.520379666564615,
                    ]
                ),
                'mosaic_id': 'P000+23',
                'related_products': 'https://vo.astron.nl/lotss_dr2/q/dlmosaic/dlmeta?ID=ivo%3A//astron.nl/%7E%3FLoTSS-DR2/P000%2B23',
                'lofar_obsids': np.array([689778]),
                'data_pid': '21.12136/ede3eff0-266a-465a-af68-6f292391485f',
                'adler32': '8eefe6e9',
            }
        ]
    elif 'P124+62' in query_string:
        result = [
            {
                'accref': 'https://vo.astron.nl/getproduct/LoTSS-DR2/P124%2B62',
                'owner': '',
                'embargo': '',
                'mime': 'image/fits',
                'accsize': np.int64(686577600),
                'centeralpha': 124.678,
                'centerdelta': 61.9895,
                'imageTitle': 'P124+62_mosaic-blanked.fits',
                'instid': 'LOFAR.HBA',
                'dateobs': 58350.241377315,
                'nAxes': 2,
                'pixelsize': np.array([9264, 9264]),
                'pixelScale': np.array([0.0004166669968981296, 0.0004166669968981296]),
                'refframe': 'ICRS',
                'wcs_equinox': 2000.0,
                'wcs_projection': 'SIN',
                'wcs_refpixel': np.array([4632.0, 4632.0]),
                'wcs_refvalues': np.array([124.678, 61.9895]),
                'wcs_cdmatrix': np.array([-0.0004166669968981296, 0.0, 0.0, 0.0004166669968981296]),
                'bandpassid': '120-168MHz',
                'bandpassunit': 'm',
                'bandpassrefval': 2.08189,
                'bandpasshi': 2.49827,
                'bandpasslo': 1.78448,
                'pixflags': 'Z',
                'coverage': np.array(
                    [
                        128.5403335917227,
                        60.00207770393527,
                        129.06122140208868,
                        63.85463521484545,
                        120.2938346473764,
                        63.85460706549431,
                        120.81483402706476,
                        60.002052894353966,
                    ]
                ),
                'mosaic_id': 'P124+62',
                'related_products': 'https://vo.astron.nl/lotss_dr2/q/dlmosaic/dlmeta?ID=ivo%3A//astron.nl/%7E%3FLoTSS-DR2/P124%2B62',
                'lofar_obsids': np.array([664568]),
                'data_pid': '21.12136/462dc630-a03c-40e9-b902-1f61b2522f42',
                'adler32': '865c7c66',
            }
        ]
    elif 'P005+21' in query_string:
        result = [
            {
                'accref': 'https://vo.astron.nl/getproduct/LoTSS-DR2/P005%2B21',
                'owner': '',
                'embargo': '',
                'mime': 'image/fits',
                'accsize': 740666880,
                'centeralpha': 5.03954,
                'centerdelta': 20.9288,
                'imageTitle': 'P005+21_mosaic-blanked.fits',
                'instid': 'LOFAR.HBA',
                'dateobs': 58668.1213425924,
                'nAxes': 2,
                'pixelsize': [9622, 9622],
                'pixelScale': [0.0004166669968981296, 0.0004166669968981296],
                'refframe': 'ICRS',
                'wcs_equinox': 2000.0,
                'wcs_projection': 'SIN',
                'wcs_refpixel': [4811.0, 4811.0],
                'wcs_refvalues': [5.03954, 20.9288],
                'wcs_cdmatrix': [-0.0004166669968981296, 0.0, 0.0, 0.0004166669968981296],
                'bandpassid': '120-168MHz',
                'bandpassunit': 'm',
                'bandpassrefval': 2.08189,
                'bandpasshi': 2.49827,
                'bandpasslo': 1.78448,
                'pixflags': 'Z',
                'coverage': [
                    7.158540286752091,
                    18.91097641889946,
                    7.2160261134913455,
                    22.92018514662428,
                    2.8626012677244765,
                    22.920179486636123,
                    2.9200990414869863,
                    18.91097090833646,
                ],
                'mosaic_id': 'P005+21',
                'related_products': 'https://vo.astron.nl/lotss_dr2/q/dlmosaic/dlmeta?ID=ivo%3A//astron.nl/%7E%3FLoTSS-DR2/P005%2B21',
                'lofar_obsids': [727372, 734329],
                'data_pid': '21.12136/85b13e46-acb7-4b58-af65-5f5e12aef98d',
                'adler32': 'dddca9d7',
            }
        ]
    elif 'P002+18' in query_string:
        result = [{
            'accref': 'https://vo.astron.nl/getproduct/LoTSS-DR2/P002%2B18',
            'owner': '',
            'embargo': '',
            'mime': 'image/fits',
            'accsize': 833955840,
            'centeralpha': 2.1779,
            'centerdelta': 18.3938,
            'imageTitle': 'P002+18_mosaic-blanked.fits',
            'instid': 'LOFAR.HBA',
            'dateobs': 58830.5910648149,
            'nAxes': 2,
            'pixelsize': [10210, 10210],
            'pixelScale': [0.0004166669968981296, 0.0004166669968981296],
            'refframe': 'ICRS',
            'wcs_equinox': 2000.0,
            'wcs_projection': 'SIN',
            'wcs_refpixel': [5105.0, 5105.0],
            'wcs_refvalues': [2.1779, 18.3938],
            'wcs_cdmatrix': [-0.0004166669968981296, 0.0, 0.0, 0.0004166669968981296],
            'bandpassid': '120-168MHz',
            'bandpassunit': 'm',
            'bandpassrefval': 2.08189,
            'bandpasshi': 2.49827,
            'bandpasslo': 1.78448,
            'pixflags': 'Z',
            'coverage': [
                4.393656223264431,
                16.25365835830698,
                4.449063431589719,
                20.508061174463997,
                359.90629143539604,
                20.50805595644316,
                359.96170949604596,
                16.253653267521848,
            ],
            'mosaic_id': 'P002+18',
            'related_products': 'https://vo.astron.nl/lotss_dr2/q/dlmosaic/dlmeta?ID=ivo%3A//astron.nl/%7E%3FLoTSS-DR2/P002%2B18',
            'lofar_obsids': [762095],
            'data_pid': '21.12136/4d6cb08b-6897-4cc9-b6b1-456690e0a141',
            'adler32': 'b56b94c3',
        }]

    return result


def _get_vo_mock(_):
    content = """<?xml-stylesheet href='/static/xsl/datalink-to-html.xsl' type='text/xsl'?>
<VOTABLE version="1.3" xmlns="http://www.ivoa.net/xml/VOTable/v1.3" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.ivoa.net/xml/VOTable/v1.3 http://vo.ari.uni-heidelberg.de/docs/schemata/VOTable-1.3.xsd">
<RESOURCE type="results">
<TABLE name="dlresponse"><DESCRIPTION>Data links for data sets.</DESCRIPTION>
<FIELD ID="ID" arraysize="*" datatype="char" name="ID" ucd="meta.id;meta.main"><DESCRIPTION>Publisher data set id; this is an identifier for the dataset in question and can be used to retrieve the data.</DESCRIPTION></FIELD>
<FIELD ID="access_url" arraysize="*" datatype="char" name="access_url" ucd="meta.ref.url"><DESCRIPTION>URL to retrieve the data or access the service.</DESCRIPTION></FIELD>
<FIELD ID="service_def" arraysize="*" datatype="char" name="service_def" ucd="meta.code"><DESCRIPTION>Identifier for the type of service if accessURL refers to a service.</DESCRIPTION></FIELD>
<FIELD ID="error_message" arraysize="*" datatype="char" name="error_message" ucd="meta.code.error"><DESCRIPTION>If accessURL is empty, this column gives the reason why.</DESCRIPTION></FIELD>
<FIELD ID="description" arraysize="*" datatype="char" name="description" ucd="meta.note"><DESCRIPTION>More information on this link</DESCRIPTION></FIELD>
<FIELD ID="semantics" arraysize="*" datatype="char" name="semantics" ucd="meta.code"><DESCRIPTION>What kind of data is linked here? Standard identifiers here include science, calibration, preview, info, auxiliary</DESCRIPTION></FIELD>
<FIELD ID="content_type" arraysize="*" datatype="char" name="content_type" ucd="meta.code.mime"><DESCRIPTION>Media type for the data returned.</DESCRIPTION></FIELD>
<FIELD ID="content_length" datatype="long" name="content_length" ucd="phys.size;meta.file" unit="byte"><DESCRIPTION>Size of the resource at access_url</DESCRIPTION><VALUES null="-1"></VALUES></FIELD>
<DATA><TABLEDATA>
<TR><TD>ivo://astron.nl/~?LoTSS-DR2/P000+23</TD><TD>https://vo.astron.nl/getproduct/LoTSS-DR2/P000+23</TD><TD></TD><TD></TD><TD>6" resolution Stokes I 120-168MHz mosaic image (fits format)</TD><TD>#this</TD><TD>image/fits</TD><TD>639106560</TD></TR>
<TR><TD>ivo://astron.nl/~?LoTSS-DR2/P000+23</TD><TD>https://vo.astron.nl/getproduct/LoTSS-DR2/P000+23?preview=true</TD><TD></TD><TD></TD><TD>Preview of the observed area</TD><TD>#preview</TD><TD>image/fits</TD><TD>120000</TD></TR>
<TR><TD>ivo://astron.nl/~?LoTSS-DR2/P000+23</TD><TD>https://lofar-webdav.grid.surfsara.nl:2881/P000+23/low-mosaic-blanked.fits</TD><TD></TD><TD></TD><TD>20" resolution Stokes I 120-168MHz mosaic image (fits format)</TD><TD>#coderived</TD><TD>image/fits</TD><TD>75000000</TD></TR>
<TR><TD>ivo://astron.nl/~?LoTSS-DR2/P000+23</TD><TD>https://lofar-webdav.grid.surfsara.nl:2881/P000+23/mosaic-weights.fits</TD><TD></TD><TD></TD><TD>6" resolution Stokes I 120-168MHz mosaic weight image (fits format)</TD><TD>#weight</TD><TD>image/fits</TD><TD>675000000</TD></TR>
<TR><TD>ivo://astron.nl/~?LoTSS-DR2/P000+23</TD><TD>https://lofar-webdav.grid.surfsara.nl:2881/P000+23/low-mosaic-weights.fits</TD><TD></TD><TD></TD><TD>20" resolution Stokes I 120-168MHz mosaic weight image (fits format)</TD><TD>#weight</TD><TD>image/fits</TD><TD>75000000</TD></TR>
<TR><TD>ivo://astron.nl/~?LoTSS-DR2/P000+23</TD><TD>https://lofar-webdav.grid.surfsara.nl:2881/P000+23/mosaic-rms.fits</TD><TD></TD><TD></TD><TD>6" resolution Stokes I 120-168MHz mosaic PYBDSF noise image (fits format)</TD><TD>#noise</TD><TD>image/fits</TD><TD>350000000</TD></TR>
<TR><TD>ivo://astron.nl/~?LoTSS-DR2/P000+23</TD><TD>https://lofar-webdav.grid.surfsara.nl:2881/P000+23/mosaic.pybdsmmask.fits</TD><TD></TD><TD></TD><TD>6" resolution Stokes I 120-168MHz mosaic PYBDSF mask image (fits format)</TD><TD>#auxiliary</TD><TD>image/fits</TD><TD>350000000</TD></TR>
<TR><TD>ivo://astron.nl/~?LoTSS-DR2/P000+23</TD><TD>https://lofar-webdav.grid.surfsara.nl:2881/P000+23/mosaic.resid.fits</TD><TD></TD><TD></TD><TD>6" resolution Stokes I 120-168MHz mosaic PYBDSF residual image (fits format)</TD><TD>#auxiliary</TD><TD>image/fits</TD><TD>675000000</TD></TR>
<TR><TD>ivo://astron.nl/~?LoTSS-DR2/P000+23</TD><TD>https://lofar-webdav.grid.surfsara.nl:2881/P000+23/fits_headers.tar</TD><TD></TD><TD></TD><TD>Headers of all the fits files (tarred).</TD><TD>#detached-header</TD><TD>application/tar</TD><TD>30720</TD></TR>
<TR><TD>ivo://astron.nl/~?LoTSS-DR2/P000+23</TD><TD>https://repository.surfsara.nl/datasets/lotss-dr2/P000-23?share-token=02786f4a-5029-b67b-424f-f7b1b4efea7b</TD><TD></TD><TD></TD><TD>Stokes I central pointing images (tarred fits format; on tape) [images.tar]</TD><TD>#progenitor</TD><TD>text/html</TD><TD>3000000</TD></TR>
<TR><TD>ivo://astron.nl/~?LoTSS-DR2/P000+23</TD><TD>https://repository.surfsara.nl/datasets/lotss-dr2/P000-23?share-token=02786f4a-5029-b67b-424f-f7b1b4efea7b</TD><TD></TD><TD></TD><TD>Central pointing direction independently calibrated measurement sets and direction dependent calibration solutions (tarred MS format and numpy arrays; on tape) [uv.tar]</TD><TD>#progenitor</TD><TD>text/html</TD><TD>3000000</TD></TR>
<TR><TD>ivo://astron.nl/~?LoTSS-DR2/P000+23</TD><TD>https://repository.surfsara.nl/datasets/lotss-dr2/P000-23?share-token=02786f4a-5029-b67b-424f-f7b1b4efea7b</TD><TD></TD><TD></TD><TD>Central pointing Stokes Q and U 20" resolution 120-168MHz image cubes with 97.6kHz planes (tarred and compressed .fits images; on tape) [stokes_A.tar]</TD><TD>#coderived</TD><TD>text/html</TD><TD>3000000</TD></TR>
<TR><TD>ivo://astron.nl/~?LoTSS-DR2/P000+23</TD><TD>https://repository.surfsara.nl/datasets/lotss-dr2/P000-23?share-token=02786f4a-5029-b67b-424f-f7b1b4efea7b</TD><TD></TD><TD></TD><TD>Central pointing 20" resolution Stokes V 120-168MHz images and Stokes Q and U 4' resolution 120-168MHz image cubes with 97.6kHz planes (tarred and compressed .fits images; on tape) [stokes_B.tar]</TD><TD>#coderived</TD><TD>text/html</TD><TD>3000000</TD></TR>
<TR><TD>ivo://astron.nl/~?LoTSS-DR2/P000+23</TD><TD>https://repository.surfsara.nl/datasets/lotss-dr2/P000-23?share-token=02786f4a-5029-b67b-424f-f7b1b4efea7b</TD><TD></TD><TD></TD><TD>Central pointing DDF-pipeline logs (tarred; on tape) [misc.tar]</TD><TD>#auxiliary</TD><TD>text/html</TD><TD>3000000</TD></TR>
<TR><TD>ivo://astron.nl/~?LoTSS-DR2/P000+23</TD><TD>https://lta.lofar.eu/Lofar?project=ALL&amp;product=all_observation_pipeline&amp;mode=query_result_page_user&amp;ObservationId=689778</TD><TD></TD><TD></TD><TD>Raw data in the Lofar Long-Term Archive</TD><TD>#progenitor</TD><TD>text/html</TD><TD>20000</TD></TR></TABLEDATA></DATA></TABLE></RESOURCE></VOTABLE>"""
    result = type('response', (), {})()
    result.close = lambda: None
    result.text = content
    return result

# https://lta.lofar.eu/Lofar?project=ALL&mode=query_result_page&product=CorrelatedDataProduct&source_pipeline=99977D292AE15486E053164A17ACDB1B&object_type=AveragingPipeline


def _map_db_query(file_name):
    results = []
    from collections import namedtuple
    CorrelatedDataProduct = namedtuple(
        'CorrelatedDataProduct',
        'data_product_id central_frequency channel_width start_time end_time ra dec file_name creator integration_interval duration file_format release_date project_name content_type content_checksum'
    )
    MinimalArtifact = namedtuple(
        'MinimalArtifact',
        'content_checksum',
    )
    with open(file_name) as f:
        content = f.read()

    from bs4 import BeautifulSoup
    from astropy.io import ascii
    soup = BeautifulSoup(content, features='lxml')
    table_header = soup.find_all('thead')
    table_body = soup.find_all('tbody')
    print(len(table_body))
    table_body_length = len(table_body)
    print(table_body[table_body_length - 1])
    file_infos = {}
    for bodies in table_body[(table_body_length - 1):]:
        file_rows = bodies.find_all('tr')
        file_name = None
        content_checksum = None

        for file_row in file_rows:
            cell_content = file_row.find_all('td')
            if cell_content[0].text == 'Filename':
                file_name = cell_content[1].text
            elif cell_content[0].text == 'Hash Md5':
                content_checksum = f'md5:{cell_content[1].text}'

        if file_name:
            file_infos[file_name] = MinimalArtifact(content_checksum)

    rows = table_body[0].findAll('tr')
    columns = table_header[0].text.split('\n')
    for index, column in enumerate(columns):
        for row in rows:

            cells = row.findAll('td')
            # print(len(cells))
            # print(cells[0])
            data_product_id = cells[6].text
            central_frequency = cells[10].text
            channel_width = cells[11].text
            start_time = cells[14].text
            end_time = cells[16].text
            ra = cells[8].text
            dec = cells[9].text
            file_name = cells[27].text
            creator = cells[3].text
            integration_interval = cells[13].text
            duration = cells[15].text
            file_format = cells[26].text
            release_date = cells[5].text
            project_name = cells[2].text
            content_type = cells[26].text

            file_info = file_infos.get(file_name)
            if file_info:
                content_checksum = file_info.content_checksum

            row_content = CorrelatedDataProduct(
                data_product_id,
                central_frequency,
                channel_width,
                start_time,
                end_time,
                ra,
                dec,
                file_name,
                creator,
                integration_interval,
                duration,
                file_format,
                release_date,
                project_name,
                content_type,
                content_checksum,
            )
            results.append(row_content)
    # print(len(rows))
    # for row in rows:
    #     cells = row.findAll('td')
    #     print(len(cells))
        # print(f'::: {row}')
    return results


def _get_db_query_mock_P164():
    return _map_db_query('/usr/src/app/lotss2caom2/lotss2caom2/tests/data/provenance/correlated_664568.htm')


def _get_db_query_mock_P005A():
    # https://lta.lofar.eu/Lofar?project=ALL&mode=query_result_page&product=CorrelatedDataProduct&pipeline_object_id=90C703E9A9705084E053144A17AC6B3C
    return _map_db_query('/usr/src/app/lotss2caom2/lotss2caom2/tests/data/provenance/correlated_734329.htm')


def _get_db_query_mock_P005B():
    # https://lta.lofar.eu/Lofar?project=ALL&mode=query_result_page&product=CorrelatedDataProduct&pipeline_object_id=8CEF25CF085D6A53E053164A17AC7525
    return _map_db_query('/usr/src/app/lotss2caom2/lotss2caom2/tests/data/provenance/correlated_727372.htm')


def _observation(collection, observation_id):
    return SimpleObservation(
        algorithm=Algorithm(name='exposure'),
        collection=collection,
        observation_id=observation_id,
    )

# if __name__ == '__main__':
   # t =  _get_db_query_mock()
   # print(t)
   # print(t[0].file_format)
