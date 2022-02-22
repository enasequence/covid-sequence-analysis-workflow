#!/usr/bin/env python3

import cx_Oracle
import os
import sys
import xmltodict


def get_sample_xml(run_id):
    client_lib_dir = os.getenv('ORACLE_CLIENT_LIB')
    if not client_lib_dir or not os.path.isdir(client_lib_dir):
        sys.stderr.write("ERROR: Environment variable $ORACLE_CLIENT_LIB must point at a valid directory\n")
        exit(1)
    cx_Oracle.init_oracle_client(lib_dir=client_lib_dir)

    era_conn = setup_connection("ora-vm-069.ebi.ac.uk", 1541, "ERAREAD", 'era_reader', 'reader')
    cursor = era_conn.cursor()

    sql = f"""
    SELECT (s.sample_xml).getClobVal() from RUN r
        join RUN_SAMPLE rs on r.run_id=rs.run_id
        join SAMPLE s on rs.sample_id=s.sample_id
    where r.run_id = '{run_id}'
    """
    cursor.execute(sql)
    result = cursor.fetchone()
    return result[0].read()


def setup_connection(host, port, service_name, user, pwd):
    connection = None
    try:
        dsn = cx_Oracle.makedsn(host, port, service_name=service_name)
        connection = cx_Oracle.connect(user, pwd, dsn, encoding="UTF-8")
        return connection
    except cx_Oracle.Error as error:
        sys.stderr.write("Could not connect to {}...\n{}\n".format(service_name, error))
        sys.exit(1)

    return connection


def get_attr(attrs, attr_name):
    try:
        return attrs[attr_name]
    except KeyError:
        return 'not provided'


def get_sample_attributes(sample_xml):
    sample_dict = xmltodict.parse(sample_xml)
    sample_acc = sample_dict['SAMPLE_SET']['SAMPLE']['@accession']

    sample_attributes_raw = {}
    if type(sample_dict['SAMPLE_SET']['SAMPLE']['SAMPLE_ATTRIBUTES']['SAMPLE_ATTRIBUTE']) is list:
        sample_attributes_raw = {attr['TAG']: attr['VALUE'] for attr in
                                 sample_dict['SAMPLE_SET']['SAMPLE']['SAMPLE_ATTRIBUTES']['SAMPLE_ATTRIBUTE']}
    else:
        attr = sample_dict['SAMPLE_SET']['SAMPLE']['SAMPLE_ATTRIBUTES']['SAMPLE_ATTRIBUTE']
        sample_attributes_raw[attr['TAG']] = attr['VALUE']

    # normalise casing and spacing of metadata fields
    sample_attributes = {k.lower().replace(' ', '_'): v for k, v in sample_attributes_raw.items()}

    attributes = {}
    attributes['isolate'] = get_attr(sample_attributes, 'isolate')
    attributes['collection_date'] = get_attr(sample_attributes, 'collection_date')

    if sample_acc[:3] == 'SRS':
        attributes['location'] = get_attr(sample_attributes, 'geo_loc_name')
    else:
        attributes['location'] = get_attr(sample_attributes, 'geographic_location_(country_and/or_sea)')
        try:
            region = sample_attributes['geographic_location_(region_and_locality)']
            if not region[:4] == 'not ':
                attributes['location'] += f":{region}"
        except KeyError:
            pass

    return attributes


def replace_header(run_id):
    sample_xml = get_sample_xml(run_id)
    sample_attributes = get_sample_attributes(sample_xml)

    # format new header
    new_header = ">CDP:{} {}".format(run_id, '|'.join([
        sample_attributes['isolate'], sample_attributes['location'],
        sample_attributes['collection_date'],
        "SARS-CoV-2 complete genome produced by the VEO task force from public INSDC data",
        "covid_sequence_analysis_workflow"
    ]))
    return new_header


with open(sys.argv[1], 'r') as fasta:
    for line in fasta:
        line = line.strip()
        try:
            if line[0] == '>':
                print(replace_header(line[1:]))
            else:
                print(line)
        except IndexError:
            pass
